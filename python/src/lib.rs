use numpy::{PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2, PyArrayMethods, IntoPyArray};
use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use sigring::{SigRing, LogSigPrep};

enum Eng { I64(SigRing<i64>), I128(SigRing<i128>) }

#[pyclass]
struct WindowSig {
    eng: Eng, prep: LogSigPrep, d: usize, w: usize, scale: f64,
    ring: Vec<i64>, head: usize, filled: usize,
}
impl WindowSig {
    fn step(&mut self, inc: &[i64]) {
        let d = self.d; let h = self.head;
        let old: Vec<i64> = self.ring[h * d..(h + 1) * d].to_vec();
        self.ring[h * d..(h + 1) * d].copy_from_slice(inc);
        self.head = (h + 1) % self.w; if self.filled < self.w { self.filled += 1; }
        match &mut self.eng { Eng::I64(e) => e.ingest(inc, &old), Eng::I128(e) => e.ingest(inc, &old) }
    }
    fn sig_vec(&self) -> Vec<f64> { match &self.eng { Eng::I64(e) => e.signature(self.scale), Eng::I128(e) => e.signature(self.scale) } }
    fn readout(&self, kind: u8) -> Vec<f64> {
        let s = self.sig_vec();
        if kind == 0 { return s; }
        let ex = sigring::tensor_log(self.d, self.n(), &s);
        match kind { 1 => ex, 2 => self.prep.lyndon_word_coefs(&ex), _ => self.prep.bracket_coords(&ex) }
    }
    fn n(&self) -> usize { match &self.eng { Eng::I64(e) => e.n, Eng::I128(e) => e.n } }
    fn out_len(&self, kind: u8) -> usize { let (d, n) = (self.d, self.n()); match kind { 0 | 1 => { let mut s = 0; let mut p = 1; for _ in 0..n { p *= d; s += p; } s } _ => self.prep.len() } }
}
fn kind_of(s: &str) -> PyResult<u8> { Ok(match s { "sig" => 0, "logsig_expanded" => 1, "logsig_words" => 2, "logsig_brackets" => 3, _ => return Err(PyValueError::new_err("readout must be sig|logsig_expanded|logsig_words|logsig_brackets")) }) }

#[pymethods]
impl WindowSig {
    /// d: dimension, n: depth, window: number of increments in the window, scale: increment integer units per 1.0 real unit, bits: 64 or 128
    #[new]
    #[pyo3(signature=(d, n, window, scale=1.0, bits=64))]
    fn new(d: usize, n: usize, window: usize, scale: f64, bits: u32) -> PyResult<Self> {
        let eng = match bits { 64 => Eng::I64(SigRing::new(d, n)), 128 => Eng::I128(SigRing::new(d, n)), _ => return Err(PyValueError::new_err("bits must be 64 or 128")) };
        Ok(Self { eng, prep: LogSigPrep::new(d, n), d, w: window, scale, ring: vec![0; d * window], head: 0, filled: 0 })
    }
    /// Push one integer increment (length d).
    fn push(&mut self, inc: PyReadonlyArray1<i64>) -> PyResult<()> {
        let s = inc.as_slice()?; if s.len() != self.d { return Err(PyValueError::new_err("bad length")); }
        self.step(s); Ok(())
    }
    /// Push a (T,d) int64 array; returns an array of readouts, one every `stride` ticks once the window is full.
    #[pyo3(signature=(incs, readout="sig", stride=1))]
    fn push_batch<'py>(&mut self, py: Python<'py>, incs: PyReadonlyArray2<i64>, readout: &str, stride: usize) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let kind = kind_of(readout)?; let a = incs.as_array(); if a.shape()[1] != self.d { return Err(PyValueError::new_err("bad width")); }
        let stride = stride.max(1); let l = self.out_len(kind);
        let mut out: Vec<f64> = Vec::new(); let mut rows = 0usize; let mut cnt = 0usize;
        let mut tmp = vec![0i64; self.d];
        py.allow_threads(|| {});
        for r in a.rows() {
            for (i, v) in r.iter().enumerate() { tmp[i] = *v; }
            self.step(&tmp);
            if self.filled == self.w { if cnt % stride == 0 { out.extend(self.readout(kind)); rows += 1; } cnt += 1; }
        }
        let arr = numpy::ndarray::Array2::from_shape_vec((rows, l), out).map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(arr.into_pyarray(py))
    }
    /// Push only (no readout) a (T,d) int64 array: pure-ingest throughput path.
    fn ingest_batch(&mut self, incs: PyReadonlyArray2<i64>) -> PyResult<()> {
        let a = incs.as_array(); let mut tmp = vec![0i64; self.d];
        for r in a.rows() { for (i, v) in r.iter().enumerate() { tmp[i] = *v; } self.step(&tmp); } Ok(())
    }
    #[pyo3(signature=(readout="sig"))]
    fn read<'py>(&self, py: Python<'py>, readout: &str) -> PyResult<Bound<'py, PyArray1<f64>>> { Ok(self.readout(kind_of(readout)?).into_pyarray(py)) }
    fn guard_ok(&self) -> bool { match &self.eng { Eng::I64(e) => e.guard_ok(), Eng::I128(e) => e.guard_ok() } }
    fn full(&self) -> bool { self.filled == self.w }
    /// Raw integer state as Python ints (for bit-exactness checks)
    fn state(&self) -> Vec<i128> { match &self.eng { Eng::I64(e) => e.t.iter().map(|&x| x as i128).collect(), Eng::I128(e) => e.t.clone() } }
}
#[pymodule]
fn _sigring(m: &Bound<'_, PyModule>) -> PyResult<()> { m.add_class::<WindowSig>()?; Ok(()) }

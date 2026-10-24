const fs=require('fs');
(async()=>{
 const {instance}=await WebAssembly.instantiate(fs.readFileSync('det.wasm'),{});
 const e=instance.exports;
 const cfgs=[[2,4,64],[3,3,32],[5,3,16],[2,5,40]];
 for(const [d,n,w] of cfgs){
   const cnt=e.run(d,n,w,200000,50000);
   const mem=new BigUint64Array(e.memory.buffer,e.outptr(),cnt);
   for(let i=0;i<cnt;i++) console.log(`${d} ${n} ${w} 64 ${(i+1)*50000} ${mem[i].toString(16).padStart(16,'0')}`);
 }
})();

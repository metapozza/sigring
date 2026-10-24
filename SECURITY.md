# Security policy

`sigring` is a numerical library with no network or file access. The relevant risks are incorrect results (for example an arithmetic or guard bug that certifies an inexact state) and memory-safety problems in the Python bindings.

Please report suspected issues privately to **metapozza@proton.me** (subject: `sigring security`) rather than in a public issue. Include the version, the inputs, and the expected versus observed output. You will get an acknowledgement within a few days. Confirmed correctness bugs are fixed in a patch release and listed in `CHANGELOG.md`.

Supported versions: the latest 0.1.x release.

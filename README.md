cif-tools
=========

The cif-tools suite of programs are tools you can use to examine
and manipulate mmCIF and PDB files.

Requirements
------------

The tools are based on [libcif++](https://github.com/PDB-REDO/libcifpp)
and the code is written in modern C++ so you need a compiler capable
of handling C++17 code.

Building
--------

Make sure you install libcif++ first before building.

After that, building should be as easy as typing:

```
./configure
make
make install
```


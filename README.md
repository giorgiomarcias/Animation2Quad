# Animation2Quad

Source code of [Animation-Aware Quadrangulation](http://vcg.isti.cnr.it/Publications/2013/MPPPS13/)

## Build

```
git clone --recursive https://github.com/giorgiomarcias/Animation2Quad.git
cd Animation2Quad/external
mkdir build_CoMISo
cd build_CoMISo
cmake ../libigl/external/CoMISo
make
```
Then just open `Animation2Quad.pro` with QtCreator and build.

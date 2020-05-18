# Install dependencies
mkdir tools
cd tools

## POA
echo -e '\n*---------------*\nInstalling POA...\n*---------------*\n'
wget 'https://ayera.dl.sourceforge.net/project/poamsa/poamsa/2.0/poaV2.tar.gz'
tar xzf poaV2.tar.gz
rm poaV2.tar.gz
cd poaV2
make poa
cd ..

## SPOA
echo -e '\n*----------------*\nInstalling SPOA...\n*----------------*\n'
git clone --recursive https://github.com/rvaser/spoa spoa
cd spoa
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -Dspoa_build_executable=ON ..
make
cd ../..


## minimap2
echo -e '\n*--------------------*\nInstalling Minimap2...\n*--------------------*\n'
git clone https://github.com/lh3/minimap2
cd minimap2 && make
cd ..

## Dense Subgraph finder
echo -e '\n*---------------*\nInstalling DSF...\n*---------------*\n'
git clone -b master --single-branch https://github.com/yana-safonova/ig_repertoire_constructor
cd ig_repertoire_constructor/
make -j 4
cd ..




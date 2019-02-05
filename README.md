# hepgen
clone the repo, then:
```
cd hepgen/HEPGenPlusPlus-1.2-RC1/build
source compile
```

To calculate the structure function at particular kinematic bin, for pi0 do:
```
./libGKPi0/gkCrossSection Q2 Xb t phi -I ../install/share/pi0_cache/tableIndex
```
or for eta:
```
./libGKPi0/gkCrossSection_eta Q2 Xb t phi -I ../install/share/eta_cache/tableIndex
```

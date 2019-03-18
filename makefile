all:
	(cd Lib && make)		
	(cd ERarefaction && make)
	(cd MetroIG && make)
	(cd MetroLogStudent && make)
	(cd MetroLogNormal && make)	
	(cd MetroSichel	&& make) 
	(cd Subsample && make)
	(cd LNAbundance && make)
	(cd IGAbundance && make)
	(cd SIAbundance && make)
	(cd LSAbundance && make)
	(cd LNRarefaction && make)
	(cd IGRarefaction && make)
	(cd SIRarefaction && make)
	(cd LSRarefaction && make)
	(cd SIShift && make)
	(cd LNShift && make)
clean:
	(cd ERarefaction && make clean)
	(cd MetroIG && make clean)
	(cd MetroLogStudent && make clean)
	(cd Lib && make clean)		
	(cd MetroLogNormal && make clean)	
	(cd MetroSichel	&& make clean) 
	(cd Subsample && make clean)
	(cd LNAbundance && make clean)
	(cd IGAbundance && make clean)
	(cd SIAbundance && make clean)
	(cd LSAbundance && make clean)
	(cd LNRarefaction && make clean)
	(cd IGRarefaction && make clean)
	(cd SIRarefaction && make clean)
	(cd LSRarefaction && make clean)
	(cd bin && rm -rf *)
	(cd SIShift && make clean)
	(cd LNShift && make clean)
install:
	cp ERarefaction/ERarefaction bin
	cp MetroIG/MetroIG bin
	cp MetroLogStudent/MetroLogStudent bin
	cp MetroLogNormal/MetroLogNormal bin	
	cp MetroSichel/MetroSichel bin
	cp Subsample/Subsample bin
	cp LNAbundance/LNAbundance bin
	cp IGAbundance/IGAbundance bin
	cp SIAbundance/SIAbundance bin
	cp LSAbundance/LSAbundance bin
	cp LNRarefaction/LNRarefaction bin
	cp IGRarefaction/IGRarefaction bin
	cp SIRarefaction/SIRarefaction bin
	cp LSRarefaction/LSRarefaction bin
	cp SIShift/SIShift bin
	cp LNShift/LNShift bin

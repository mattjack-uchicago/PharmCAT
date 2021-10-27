MAKEFLAGS = -j1

# Determine this makefile's path.
# Be sure to place this BEFORE `include` directives, if any.
THIS_FILE := $(lastword $(MAKEFILE_LIST))

ifdef PHARMCAT_DATA_DIR
  dataDir := $(PHARMCAT_DATA_DIR)
else
  dataDir := build
endif

ifeq ($(OS),Windows_NT)
	GRADLE_CMD := cmd /c gradlew.bat --console=plain
else
	GRADLE_CMD := ./gradlew --console=plain
endif


.PHONY: updateData
updateData: clean
	${GRADLE_CMD} updateData
	mv src/main/resources/org/pharmgkb/pharmcat/definition/alleles/pharmcat_positions.* .
	${GRADLE_CMD} updateExample
	mv docs/examples/pharmcat_positions.matcher.html    docs/examples/pharmcat.example.matcher.html
	mv docs/examples/pharmcat_positions.matcher.json    docs/examples/pharmcat.example.matcher.json
	mv docs/examples/pharmcat_positions.phenotyper.json docs/examples/pharmcat.example.phenotyper.json
	mv docs/examples/pharmcat_positions.report.html     docs/examples/pharmcat.example.report.html
	mv docs/examples/pharmcat_positions.report.json     docs/examples/pharmcat.example.report.json


.PHONY: docker
docker: clean
	${GRADLE_CMD} shadowJar
	mv build/libs/pharmcat-`git describe --tags | sed -r s/^v//`-all.jar build/pharmcat.jar
	docker build -t pcat .


.PHONY: scriptPkg
scriptPkg:
	rm -rf build/preprocessor
	mkdir -p build/preprocessor
	cp src/scripts/preprocessor/*.txt build/preprocessor
	cp src/scripts/preprocessor/*.py build/preprocessor
	cp pharmcat_positions.vcf* build/preprocessor
	cp PharmCAT.wiki/Preprocessing-VCF-Files-for-PharmCAT.md build/preprocessor/README.md
	cd build; tar -czvf preprocessor.tar.gz preprocessor


.PHONE: updateDataFromScratch
updateDataFromScratch: docker updateData


.PHONY: tests
tests: exactVcfTests exactVcfMissingTests


.PHONY: vcfTests
vcfTests: clean
	rm -rf ${dataDir}/testVcf ${dataDir}/autogeneratedTestResults
	src/scripts/vcf_generator/generate_vcf_test_data.sh
	${GRADLE_CMD} testAutogeneratedVcfs
	cd ${dataDir}; zip -r vcfTest-`date +'%Y-%m-%d'`.zip autogeneratedTestResults

.PHONY: vcfMissingTests
vcfMissingTests: clean
	rm -rf ${dataDir}/testVcf ${dataDir}/autogeneratedTestResults
	src/scripts/vcf_generator/generate_vcf_test_data.sh -m
	${GRADLE_CMD} testAutogeneratedVcfs
	cd ${dataDir}; zip -r vcfTest-`date +'%Y-%m-%d'-missing`.zip autogeneratedTestResults


.PHONY: exactVcfTests
exactVcfTests: clean
	rm -rf ${dataDir}/testVcf ${dataDir}/autogeneratedTestResults
	src/scripts/vcf_generator/generate_vcf_test_data.sh
	${GRADLE_CMD} testAutogeneratedVcfsExactMatchOnly
	cd ${dataDir}; zip -r vcfTest-`date +'%Y-%m-%d'`-exact.zip autogeneratedTestResults

.PHONY: exactVcfMissingTests
exactVcfMissingTests: clean
	rm -rf ${dataDir}/testVcf ${dataDir}/autogeneratedTestResults
	src/scripts/vcf_generator/generate_vcf_test_data.sh -m
	${GRADLE_CMD} testAutogeneratedVcfsExactMatchOnly
	cd ${dataDir}; zip -r vcfTest-`date +'%Y-%m-%d'`-exact-missing.zip autogeneratedTestResults


.PHONY: fuzzyVcfTests
fuzzyVcfTests: clean
	rm -rf ${dataDir}/testVcf ${dataDir}/autogeneratedTestResults
	src/scripts/vcf_generator/generate_vcf_test_data.sh
	${GRADLE_CMD} testAutogeneratedVcfsFuzzyMatch
	cd ${dataDir}; zip -r vcfTest-`date +'%Y-%m-%d'`-fuzzy.zip autogeneratedTestResults

.PHONY: fuzzyVcfMissingTests
fuzzyVcfMissingTests: clean
	rm -rf ${dataDir}/testVcf ${dataDir}/autogeneratedTestResults
	src/scripts/vcf_generator/generate_vcf_test_data.sh -m
	${GRADLE_CMD} testAutogeneratedVcfsFuzzyMatch
	cd ${dataDir}; zip -r vcfTest-`date +'%Y-%m-%d'`-fuzzy-missing.zip autogeneratedTestResults


.PHONY: clean
clean:
	${GRADLE_CMD} clean


.PHONY: release
release:
	yarn release
	@echo "Updating main branch..."
	git checkout main
	git pull
	git rebase origin/development
	git push
	# switching back to development
	git checkout development
	@echo "\nDone."


.PHONY: dockerRelease
.ONESHELL:
dockerRelease: docker
	version=`git describe --tags | sed -r s/^v//`
	docker tag pcat pgkb/pharmcat:$${version}
	docker push pgkb/pharmcat:$${version}
	docker tag pcat pgkb/pharmcat:latest
	docker push pgkb/pharmcat:latest

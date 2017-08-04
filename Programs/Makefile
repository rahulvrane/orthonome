#**********************************************#
# Installation script adopted from MSOAR2 installation

CC = g++

all: mapToPhylip familySeparator addOutgroup getTopHits postprocessing removeInparalogs normalizeScores tagGenerator treeTran

mapToPhylip:
	$(CC)  MapToPhylip.cpp -o mapToPhylip

familySeparator:
	$(CC)  FamilyAssembler.cpp -o familyAssembler

addOutgroup:
	$(CC)  AddOutgroup.cpp -o addOutgroup

treeTran:
	$(CC)  TreeTran.cpp -o treeTran

postprocessing:
	$(CC)  Postprocessing.cpp -o postprocessing

removeInparalogs:
	$(CC)  RemoveInparalogs2.cpp -o removeInparalogs2

normalizeScores:
	$(CC)  NormalizeScores.cpp -o normalizeScores

tagGenerator:
	$(CC)  TagGenerator.cpp -o tagGenerator

clean: 
	rm -f *~ mapToPhylip familyAssembler addOutgroup getTopHits postprocessing removeInparalogs2 normalizeScores tagGenerator treeTran

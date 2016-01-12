#include <iostream>
#include <vector>
#include "simplicialChainComplex.h"
using namespace std;

int main (int argc, char *argv[]) {

	// read	
	simplicialPolyhedron P;
	P.read(cin);
	simplicialChainComplex S(P);

	// prepare
	int * deletedSimplexes = (int *) malloc(0*sizeof(int));
	int numDeletedSimplexes = 0;

	// cell complex construction
	for (int i = S.dim(); i > 0; i--) {
cerr << i << endl;
		sparseMatrix L,D,U,rP,cP;
		sparseMatrix d = S.boundaryOperator(i).transpose();
d.print_octave(cerr);
for (int l = 0; l < numDeletedSimplexes; l++)
	cerr << " " << deletedSimplexes[l];
cerr << endl;
		d.deleteRows(numDeletedSimplexes,deletedSimplexes);
d.print_octave(cerr);
		d.LDU_efficient(L,D,U,rP,cP);
		
		// interprete
		numDeletedSimplexes = 0;
		deletedSimplexes = (int *) realloc(deletedSimplexes, L.size(2)*sizeof(int));
		sparseMatrix Lt = L.transpose();
		sparseMatrix rPt = rP.transpose();
		for (int col = 0; col < Lt.size(1); col++) {
			int numValuesInCol = Lt.numValuesInRow(col);
			if (numValuesInCol == 1) {
				// cell retraction
				int cell1 = rPt.getCols()[Lt.getCols()[Lt.getRows()[col]]];
				int face = cP.getCols()[col];
cerr << "Cells " << cell1 << " retracted throught " << face << endl;

				deletedSimplexes[numDeletedSimplexes] = face;
				numDeletedSimplexes++;

			} else if (numValuesInCol == 2) {
				// cell joining
				int cell1 = rPt.getCols()[Lt.getCols()[Lt.getRows()[col]]];
				int cell2 = rPt.getCols()[Lt.getCols()[Lt.getRows()[col]+1]];
				int face = cP.getCols()[col];

cerr << "Join between cells " << cell1 << " and " << cell2 << " throught " << face << endl;
				deletedSimplexes[numDeletedSimplexes] = face;
				numDeletedSimplexes++;

			} else {
				// cell spliting
			}

		// arrange
		}
	}

}
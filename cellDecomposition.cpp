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
cerr << "Simplexes of dimension " << i << endl;
		sparseMatrix L,D,U,rP,cP;
		sparseMatrix d = S.boundaryOperator(i).transpose();
		int numAliveCells = (d.size(1) - numDeletedSimplexes);
		int * aliveCellsIndex = (int *) malloc(numAliveCells*sizeof(int));

		for (int j = 0, ac = 0 , dc = 0; j < d.size(1) ; j++)
			if ( dc < numDeletedSimplexes && j == deletedSimplexes[dc]) {
				dc++;
			} else {
				aliveCellsIndex[ac] = j;
				ac++;
			}

		d.deleteRows(numDeletedSimplexes,deletedSimplexes);
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
				int cell1 = aliveCellsIndex[rPt.getCols()[Lt.getCols()[Lt.getRows()[col]]]];
				int face = cP.getCols()[col];
cerr << "Cells " << cell1 << " retracted throught " << face << endl;

				deletedSimplexes[numDeletedSimplexes] = face;
				numDeletedSimplexes++;

			} else if (numValuesInCol == 2) {
				// cell joining
				int cell1 = aliveCellsIndex[rPt.getCols()[Lt.getCols()[Lt.getRows()[col]]]];
				int cell2 = aliveCellsIndex[rPt.getCols()[Lt.getCols()[Lt.getRows()[col]+1]]];
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
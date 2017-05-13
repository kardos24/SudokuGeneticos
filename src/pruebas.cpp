#include <fstream>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <ga/GASimpleGA.h>
#include <ga/GA1DArrayGenome.h> //genoma-> array de enteros (dim. 1) alelos
#include <ga/ga.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>

using namespace std;

//Definimos una estructura para guardar los datos.
struct plantilla {
	int tam; // tamaño del sudoku
	int *fijo; // valores fijos iniciales
};

char* nombreFichero = "Casos/Sudoku-A1.txt";

float Objective(GAGenome &); //Funcion objetivo
void InicioSudoku(GAGenome& g);
int MutacionSudoku(GAGenome& g, float pmut);
int CruceSudoku(const GAGenome& p1, const GAGenome & p2, GAGenome* c1,
		GAGenome* c2);
void leerSudoku(struct plantilla *S, char *nombreF);
void mostrarSudoku(struct plantilla *S);

int main() {

	// Leemos el ficheroç
	// nombreFichero = argv[1];
	struct plantilla *sudoku = new plantilla;
	cout << "Empezar a leer el fichero " << nombreFichero << endl;
	leerSudoku(sudoku, nombreFichero);
	cout << "Terminar de leer el fichero" << endl;
	cout << "Sudoku con los valores iniciales" << endl;
	mostrarSudoku(sudoku);

	// Declaramos variables para los parametros del GA y las inicializamos
	int popsize = 60;		// atoi(argv[2]);
	int ngen = 8000;			// atoi(argv[3]);
	float pmut = 0.6;		// atoi(argv[4]);
	float pcross = 0.8;		// atoi(argv[5]);

	GAAlleleSet<int> alelos;
	for (int i = 1; i <= sudoku->tam; i++)
		alelos.add(i);

	//Creamos la estructura del GENOMA, un array de bits.
	GA1DArrayAlleleGenome<int> genome(sudoku->tam, alelos, Objective, sudoku);
//	genome.initializer(GA1DArrayAlleleGenome< vector<int> >::OrderedInitializer);
//    genome.mutator(GA1DArrayAlleleGenome< vector<int> >::SwapMutator);
//    genome.crossover(GA1DArrayAlleleGenome< vector<int> >::OrderCrossover);
	genome.initializer(InicioSudoku);
	genome.crossover(CruceSudoku);
	genome.mutator(MutacionSudoku);
//	genome.crossover(GA1DArrayAlleleGenome<int>::OnePointCrossover);

	GASimpleGA ga(genome);

	ga.minimaxi(-1);
	ga.populationSize(popsize);
	ga.nGenerations(ngen);
	ga.pMutation(pmut);
	ga.pCrossover(pcross);
	ga.evolve(1);
	cout << ga.statistics().minEver() << endl;		//Numero de desplazamientos

	return 0;
}

// Funcion objetivo
float Objective(GAGenome& g) {
	GA1DArrayAlleleGenome<int> & genome = (GA1DArrayAlleleGenome<int> &) g;

	// Obtenemos la estructura auxiliar
	struct plantilla * plantilla1;
	plantilla1 = (struct plantilla *) genome.userData();
	int numFallos = 0;
	int sq = sqrt(plantilla1->tam);

	for(int i = 0; i < plantilla1->tam; i++){	//Para cada fila, cada columna y cada cuadrante
		int auxfi[plantilla1->tam];
		int auxco[plantilla1->tam];
		int auxcu[plantilla1->tam];
		int por = i/sq;
		por *= sq;
		for(int x = 0; x < plantilla1->tam; x++){
			auxfi[x] = 0;
			auxco[x] = 0;
			auxcu[x] = 0;
		}
		for(int j = 0; j < plantilla1->tam; j++){
			auxfi[genome.gene(i*plantilla1->tam+ j)]+=1;	// Coordenadas i,j para recorrer las filas
			auxco[genome.gene(j*plantilla1->tam+ i)]+=1;	// Coordenadas i,j para recorrer las columnas
			int icua = j/sq + por;	// Coordenada i para el cuadrante i,j
			int jcua = j%sq + por;	// Coordenada j para el cuadrante i,j
			auxcu[genome.gene(icua*plantilla1->tam + jcua)]+=1;
		}
		for(int y = 0; y < plantilla1->tam; y++){
			if(auxfi[y] == 0) numFallos++;
			if(auxco[y] == 0) numFallos++;
			if(auxcu[y] == 0) numFallos++;
		}
	}

	return 0.0;
}

void InicioSudoku(GAGenome& g) {

	// Hacer el casting correspondiente para obtener genome
	GA1DArrayAlleleGenome<int> & genome = (GA1DArrayAlleleGenome<int> &) g;

	// Obtenemos la estructura auxiliar
	struct plantilla * plantilla1;
	plantilla1 = (struct plantilla *) genome.userData();

	int aux[plantilla1->tam];

	for (int f = 0; f < plantilla1->tam; f++) { // por cada fila

		// se inicializa a 0 la estructura auxiliar

		for (int j = 0; j < plantilla1->tam; j++)
			aux[j] = 0;

		// se inicializa la fila de 1 a tam sin repetidos

		for (int j = 1; j <= plantilla1->tam; j++) {
			int v = GARandomInt(0, plantilla1->tam - 1);
			while (aux[v] != 0)
				v = (v + 1) % plantilla1->tam;
			aux[v] = j;
		}

		// se colocan los fijos en su lugar

		int i = 0;

		while (i < plantilla1->tam) {

			while ((plantilla1->fijo[(f * plantilla1->tam) + i] == 0)
					&& (i < plantilla1->tam))
				i++;

			if (i < plantilla1->tam) { // se encuentra un fijo en la plantilla

				// se busca el fijo en aux

				bool encontrado = false;
				for (int j = 0; (j < plantilla1->tam) && (!encontrado); j++)
					if (aux[j] == plantilla1->fijo[(f * plantilla1->tam) + i]) {
						encontrado = true;
						aux[j] = aux[i];
					}
				// se pone el fijo en su sitio

				aux[i] = plantilla1->fijo[(f * plantilla1->tam) + i];
			}
			i++;

		}

		// se copia la fila en el genoma

		for (int c = 0; c < plantilla1->tam; c++)
			genome.gene((f * plantilla1->tam) + c, aux[c]);
	}
}

// Comprueba si una columna del sudoku tiene valores repetidos
// check contiene la cuenta de valores repetidos

bool checkColumna(int col[], int * check, int tam) {
	bool repe = false;

	for (int i = 0; i < tam; i++)
		check[i] = 0;

	for (int i = 0; i < tam; i++)
		check[col[i] - 1]++;
	for (int i = 0; i < tam; i++)
		if (check[i] > 1)
			repe = true;

	return repe;
}

// Si un gen debe mutar elegimos si lo hace en fila o columna. En fila es intercambio
// en columna es mas heurístico: cambio un valor repetido por uno que no exista en la columna

int MutacionSudoku(GAGenome& g, float pmut) {

	// Hacer el casting correspondiente para obtener genome
	GA1DArrayAlleleGenome<int> & genome = (GA1DArrayAlleleGenome<int> &) g;

	// Obtenemos la estructura auxiliar
	struct plantilla * plantilla1;
	plantilla1 = (struct plantilla *) genome.userData();
	int nmut = 0;
	int aux;
	int fil;
	bool fila;

	int caux[plantilla1->tam];
	int *checkC = new int[plantilla1->tam];

	if (pmut <= 0.0)
		return 0;

	for (int f = 0; f < plantilla1->tam; f++)
		for (int c = 0; c < plantilla1->tam; c++)
			if (plantilla1->fijo[(f * plantilla1->tam) + c] == 0) { // si no es fijo
				if (GAFlipCoin(pmut)) { // si toca mutar
					if (GAFlipCoin(0.5))
						fila = true; // cambiamos fila
					else
						fila = false; // cambiamos columna

					if (!fila) { // mutamos columna y arreglamos fila

						for (int j = 0; j < plantilla1->tam; j++)
							caux[j] = genome.gene((j * plantilla1->tam) + c); // obtenemos la columna
						if (checkColumna(caux, checkC, plantilla1->tam)) { // si hay repetidos en la columna
							int v1 = GARandomInt(0, plantilla1->tam - 1); // v1 es valor repetido
							while (checkC[v1] <= 1)
								v1 = (v1 + 1) % plantilla1->tam;
							v1++;
							int v2 = GARandomInt(0, plantilla1->tam - 1); // v2 es un valor que no existe
							while (checkC[v2] != 0)
								v2 = (v2 + 1) % plantilla1->tam;
							v2++;
							//buscamos donde está v1 y se cambia por v2

							bool encontrado = false;
							for (int j = 0; j < plantilla1->tam && !encontrado;
									j++)
								if ((plantilla1->fijo[j * (plantilla1->tam) + c]
										== 0)
										&& (genome.gene(
												j * (plantilla1->tam) + c) == v1)) {
									encontrado = true;
									genome.gene((j * plantilla1->tam) + c, v2);
									fil = j;  // v1 esta en fil
								}
							//arreglamos la fila fil donde está v1: buscamos v2 y ponemos v1

							int col = (c + 1) % plantilla1->tam;
							while (genome.gene((fil * plantilla1->tam) + col)
									!= v2)
								col = (col + 1) % plantilla1->tam;
							if (plantilla1->fijo[(fil * plantilla1->tam) + col]
									== 0) { // si v2 no es fijo en fil se lleva a cabo la mutación
								nmut++;
								genome.gene((fil * plantilla1->tam) + col, v1);
							} else { // si es fijo se deshace la mutacion
								genome.gene((fil * plantilla1->tam) + c, v1);
							}

						} // fin de si hay repetidos

					} else { // mutamos fila: cambiamos f,c por f,v1
						int v1 = (c + 1) % plantilla1->tam;
						while ((plantilla1->fijo[(f * plantilla1->tam) + v1]
								!= 0))
							v1 = (v1 + 1) % plantilla1->tam; //busco un no fijo
						aux = genome.gene((f * plantilla1->tam) + c);
						genome.gene((f * plantilla1->tam) + c,
								genome.gene((f * plantilla1->tam) + v1));
						genome.gene((f * plantilla1->tam) + v1, aux);
						nmut++;
					}
				} // si toca mutar
			} // si no fijo

	return nmut;
	return 1;
}

// Cruce por un punto que coincide con final de una fila

int CruceSudoku(const GAGenome& p1, const GAGenome & p2, GAGenome* c1,
		GAGenome* c2) {

	// Hacer el casting correspondiente para obtener m y p
	GA1DArrayAlleleGenome<int> & m = (GA1DArrayAlleleGenome<int> &) p1;
	GA1DArrayAlleleGenome<int> & p = (GA1DArrayAlleleGenome<int> &) p2;

	struct plantilla * plantilla1 = (struct plantilla *) m.userData();
	int n = 0;

	// buscamos un punto de cruce por filas

	int punto1 = GARandomInt(0, m.length());
	while ((punto1 % plantilla1->tam) != 0)
		punto1 = (punto1 + 1) % plantilla1->tam;
	int punto2 = m.length() - punto1;

	// el hijo 1 hereda la primera parte del padre 1 y la segunda parte del
	// padre 2.

	if (c1) {
		// Hacer el casting correspondiente para obtener h1
		GA1DArrayGenome<int> & h1 = (GA1DArrayGenome<int> &) c1;

		h1.copy(m, 0, 0, punto1);
		h1.copy(p, punto1, punto1, punto2);
		n++;
	}

	// el hijo 2 hereda la primera parte del padre 2 y la segunda parte del
	// padre 1.

	if (c2) {
		// Hacer el casting correspondiente para obtener h2
		GA1DArrayGenome<int> & h2 = (GA1DArrayGenome<int> &) c2;

		h2.copy(p, 0, 0, punto1);
		h2.copy(m, punto1, punto1, punto2);
		n++;
	}

	return n;

}

// Lectura del sudoku inicial en S
void leerSudoku(struct plantilla *S, char *nombreF) {

	//Leemos los Mochila de fichero.
	ifstream f(nombreF);

	//Leemos el capacidad del que disponemos
	f >> S->tam;

	//Inicilizamos el sudoku
	S->fijo = new int[S->tam * S->tam];

	//Leemos el sudoku inicial
	for (int i = 0; i < S->tam * S->tam; i++)
		f >> S->fijo[i];

	f.close();
	cout << "El tamaño del sudoku es: " << S->tam << endl;
}

// Mostrar la información del sudoku
void mostrarSudoku(struct plantilla *S) {
	for (int i = 0; i < S->tam * S->tam; i++) {
		cout << S->fijo[i];
		if ((i + 1) % S->tam == 0) {
			cout << endl;
		} else {
			cout << " ";
		}
	}
}

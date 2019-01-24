#include <mpi.h>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#define INF 1000000000
#define abs(x) ((x) < 0 ? -(x) : (x))


using namespace std;




static int N = 300;									// ilosc podzialow boku kwadratu
static double dokladnosc = 1e-5;					// dokladnosc obliczenia
static double h = 1./N;								// odleglosc oczek siatki
static double alfa = 0.5;							// wspolczynnik relaksacji



// funkcja zwracajaca czas z duza dokladnoscia:
//
double GetTickCount(void) 
{
	struct timespec now;
	if (clock_gettime(CLOCK_MONOTONIC, &now))
		return 0;
	return now.tv_sec * 1000.0 + now.tv_nsec / 1000000.0;
}


// zamiana koloru HSV na kolor RGB:
//
void hsv2rgb(double hue, double sat, double val, double &red, double &grn, double &blu) { 

	double i, f, p, q, t;

	if (val == 0) {
		red = 0;
		grn = 0;
		blu = 0;
	}
	else {
		hue /= 60;
		i = floor(hue);
		f = hue - i;
		p = val*(1 - sat);
		q = val*(1 - (sat*f));
		t = val*(1 - (sat*(1 - f)));
		if (i == 0) { red = val; grn = t; blu = p; }
		else if (i == 1) { red = q; grn = val; blu = p; }
		else if (i == 2) { red = p; grn = val; blu = t; }
		else if (i == 3) { red = p; grn = q; blu = val; }
		else if (i == 4) { red = t; grn = p; blu = val; }
		else if (i == 5) { red = val; grn = p; blu = q; }
	}
}


// kolorowy wykres danego pola:
//
void rysuj_kolorowy_wykres(string s, double **pole) {


	unsigned char *img = new unsigned char[54 + 3 * N*N];	
	
	int i, j;
	double mn = INF, mx = -INF, MN = INF;



	int w, h, x, y;
	w = N;
	h = N;

	FILE *f;

	int filesize = 54 + 3 * w*h;  //w is your image width, h is image height, both int




								  // znajdz ekstremalne wartosci pola
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			mn = min(mn, pole[i][j]),
			mx = max(mx, pole[i][j]);



	// czerwony kolor - maksymalna wartość
	// fioletowy kolor - minimalna wartość

	double r, g, b, hue;

	for (j = 0; j < N; j++)
		for (i = 0; i < N; i++) {

			hue = pole[i][j];
			hue = max(hue, mn);
			hue = min(hue, mx);

			hue = 360 - 360 * (hue - mn) / (mx - mn);
			hue *= 0.8;
			//hue = 360 - 360 * (hue - mn) / (mx - mn);

			while (hue > 360)
				hue -= 360;
			while (hue < 0)
				hue += 360;


			hsv2rgb(hue, 1, 1, r, g, b);

			r *= 255;
			g *= 255;
			b *= 255;


			x = i; y = (h - 1) - j;
			img[(x + y*w) * 3 + 2] = (unsigned char)(r);
			img[(x + y*w) * 3 + 1] = (unsigned char)(g);
			img[(x + y*w) * 3 + 0] = (unsigned char)(b);
		}





	unsigned char bmpfileheader[14] = { 'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0 };
	unsigned char bmpinfoheader[40] = { 40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0 };
	unsigned char bmppad[3] = { 0,0,0 };

	bmpfileheader[2] = (unsigned char)(filesize);
	bmpfileheader[3] = (unsigned char)(filesize >> 8);
	bmpfileheader[4] = (unsigned char)(filesize >> 16);
	bmpfileheader[5] = (unsigned char)(filesize >> 24);

	bmpinfoheader[4] = (unsigned char)(w);
	bmpinfoheader[5] = (unsigned char)(w >> 8);
	bmpinfoheader[6] = (unsigned char)(w >> 16);
	bmpinfoheader[7] = (unsigned char)(w >> 24);
	bmpinfoheader[8] = (unsigned char)(h);
	bmpinfoheader[9] = (unsigned char)(h >> 8);
	bmpinfoheader[10] = (unsigned char)(h >> 16);
	bmpinfoheader[11] = (unsigned char)(h >> 24);

	f = fopen(s.c_str(), "wb");
	fwrite(bmpfileheader, 1, 14, f);
	fwrite(bmpinfoheader, 1, 40, f);
	for (i = 0; i < h; i++)
	{
		fwrite(img + (w*(h - i - 1) * 3), 3, w, f);
		fwrite(bmppad, 1, (4 - (w * 3) % 4) % 4, f);
	}
	
	fclose(f);
	
	delete [] img;
}
	
	
	
	
	
// obliczenia na watku nr q az do zbieznosci:
//
void foo( int q, int MODS, int size, double **u) {
	// 		  ^		  ^
	// numer watku oraz ilosc watkow
	
	
	MPI_Status status;
	MPI_Request request;
	int i, j, k;
	double u_p, eps, eps_rest;
	
	
	
	
	// stworz tablice pomocnicza **u_new, zapisujaca nastepny krok iteracji, nie nadpisujac poprzedniego:
	//
	double **u_new = new double * [size+2];
	
	for( i = 0; i <= size+1; i++)
		u_new[i] = new double[N];
		

	
	
	
	
	
	// pierwsze "przyblizenie" rozwiazania **u wraz z warunkami brzegowymi:
	//
	for( i = 0; i <= size+1; i++)
		for( j = 0; j < N; j++)
			u[i][j] = 0;
			

	
	

	eps = INF;
	// wykonaj k iteracji:
	//
	for( k = 0; sqrt(eps) > dokladnosc; k++) {
		
		
		// wyslij watkom poprzedniemu i nastepnemu gorny i dolny wiersz LICZONEJ CZESCI tablicy:
		//
		
		if( q-1 >= 0)
			MPI_Isend ( u[1], N, MPI_DOUBLE, q-1, 0, MPI_COMM_WORLD, &request);
		if( q+1 < MODS)
			MPI_Isend ( u[size], N, MPI_DOUBLE, q+1, 0, MPI_COMM_WORLD, &request);
		
		
		

		// otrzymaj gorny i dolny wiersz od watkow poprzedniego i nastepnego: (TYLKO DO ODCZYTU)
		//
		
		if( q-1 >= 0)
			MPI_Recv ( u[0], N, MPI_DOUBLE, q-1, 0, MPI_COMM_WORLD, &status );
		if( q+1 < MODS)
			MPI_Recv ( u[size+1], N, MPI_DOUBLE, q+1, 0, MPI_COMM_WORLD, &status );
		
	
	
	
	

		eps = 0;
	
		// wykonaj 1 iteracje na kazdym "wewnetrznym" elemencie tablicy:
		//
		for( i = 1; i <= size; i++)
			for( j = 1; j <= N-2; j++) {
				
				u_p = ( h*h + u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1]) / 4.;
				
				u_new[i][j] = alfa * u[i][j] + (1-alfa) * u_p;
				
				eps += ( u[i][j] - u_new[i][j]) * ( u[i][j] - u_new[i][j]);
			}


		
		
		// watek glowny sumuje eps nadeslane ze wszystkich watkow:
		//
		if( q != 0) {
			
			// najpierw wyslij glownemu watkowi swoja czesc sumy:
			//
			MPI_Ssend ( &eps, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			
			// nastepnie odbierz od niego sume wszystkich:
			//
			MPI_Recv ( &eps, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		}
		else {
			
			// najpierw odbierz eps od wszystkich watkow i zsumuj:
			//
			for( i = 1; i < MODS; i++) {
				
				MPI_Recv ( &eps_rest, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
				eps += eps_rest;
			}
			
			// nastepnie odeslij kazdemu watkowi sume:
			//
			for( i = 1; i < MODS; i++)
				MPI_Send ( &eps, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}
	
	
	
	
	
		// **u_new := **u
		//
		for( i = 1; i <= size; i++)
			for( j = 1; j <= N-2; j++)
				u[i][j] = u_new[i][j];
		
		
		
		
		if( k % 1000 == 0 && k > 0) {
			
			//printf( "watek nr %d policzony w %d/10\n", q+1, k/1000);
			if( q == 0)
				printf("iteracje = %d, epsilon = %.6E\n", k, sqrt(eps));
		}
		
	}
	




	// pozbadz sie tablicy pomocniczej:
	//
	for( i = 0; i <= size+1; i++)
		delete [] u_new[i];
		
	delete [] u_new;
	




	
	// wyslij wyniki glownemu watkowi (nr 1):
	//
	if( q != 0)		// nie wysylaj do siebie ( watek 1 -/-> watek 1)
		for( i = 1; i <= size; i++) 
			MPI_Send ( u[i], N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
}






int main ( int argc, char *argv[] ) {
	
	int id;
	int MODS;
	MPI_Status status;
	MPI_Request request;
  
  
	//  Initialize MPI.
	MPI_Init ( &argc, &argv );
	
	
	//  Get the number of processes.
	//
	MPI_Comm_size ( MPI_COMM_WORLD, &MODS );
	//
	//  Get the individual process ID.
	//
	MPI_Comm_rank ( MPI_COMM_WORLD, &id );
	
	
	

	printf( "zaczynam watek nr %d\n", id+1);
	
	
	
	
	
	int i, j, k, q;
	double t;
	
	



	// tablica rozmiarow tablic, ktore otrzymaja poszczegolne watki: (znana wszystkim watkom)
	//
	int *size = new int[MODS];
	

	// ustal rozmiary fragmentow tablic, ktore zostana przekazane poszczegolnym watkom:
	//
	for( q = 0; q < MODS; q++)
		size[q] = (N-2) / MODS;
	
	// pierwsze kilka watkow dostanie o 1 wiekszy rozmiar niz pozostale:
	//
	for( q = 0; q < (N-2) - (N-2) / MODS * MODS; q++)
		size[q]++;
	
	
	
	


	// wypisz sposob rozdzielenia tablicy na watki:
	//
	if( id == 0) {
		

		printf( "sposob podzialu wierszy tablicy na watki:\n[");
		for( q = 0; q < MODS; q++)
			printf( "%d ", size[q]);
		printf("]\n");
	}
		
	// zacznij mierzyc czas wykonania programu:
	//
	if( id == 0) {

		t = GetTickCount();
	}
	
	
	
	



	

	// stworz tablice **u dla watku, powiekszona o 2 wiersze (gorny i dolny):
	//
	double **u = new double * [size[id]+2];

	for( i = 0; i <= size[id]+1; i++)
		u[i] = new double[N];
	
	
	
	// glowna czesc programu - obliczenia:
	//
	foo( id, MODS, size[id], u);





	


	// stworz tablice dynamiczna do prezentacji wynikow:
	//
	if (id == 0) {
		
		double **v;
	
		v = new double * [N];
	
		for( i = 0; i < N; i++)
			v[i] = new double[N];
	


		// zbierz wyniki ze wszystkich watkow:
		//
		k = 1;
		for( j = 0; j < size[0]; j++)	// najpierw z watku nr 1 (glownego)
			copy( u[j], u[j]+N, v[k++]);

		
		for( i = 1; i < MODS; i++)		// potem z pozostalych watkow
			for( j = 0; j < size[i]; j++)
				MPI_Recv ( v[k++], N, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status );



		rysuj_kolorowy_wykres( "laplasjan.bmp", v);
		
		
		
		for( i = 0; i < N; i++)
			delete [] v[i];
	
	
		delete [] v;
	}
	

	
	
	
	
	
	
	if( id == 0) {
		
		printf( "czas wykonania programu = %lfs\n", (GetTickCount()-t) / 1000);
		
		FILE *f = fopen( "result.txt", "a");
		
		fprintf( f, "czas wykonania programu na liczbie watkow %d wynosi %lfs\n", MODS, (GetTickCount()-t) / 1000);
		fclose(f);
	}
	
	
	
	

	
	
	for( i = 0; i <= size[id]+1; i++)
		delete [] u[i];
		
	
	delete [] u;
	delete [] size;
		
	printf( "terminuje watek nr %d\n", id+1);

	MPI_Finalize ( );
	


	return 0;
}

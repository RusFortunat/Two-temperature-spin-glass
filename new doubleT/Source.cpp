#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <time.h>
#include <random>
#include <cmath>
#include <fstream>
#include <sstream>

using namespace std;

// ==================================================
// Program description
// ==================================================
//
//	This is a simulation of the KLS System, part of which hold at critical temperature and other part at high temperature.
//	Specifically, I model spin exchange dynamics in the lattice gas with one spin per one lattice site. At the beginning 
//	+1/-1 spins are distributed randomly in the lattice. After switching on the infinite drive field and Ising interactions,
//	system relaxes to its steady state. The details of relaxation could be studied through two-times correlation function 
//	and computing it is the purpose of this program.

int main() {

	// ==================================================
	// Simulation parameters
	// ==================================================

	int i, j, n, dt, delta, E1, iPrev, iPPrev, iNext, iNNext, jPrev, jNext, dice;
	double r, w, W;
	int Lx = 40;  // Width
	int Ly = (Lx*Lx*Lx);  // Length
	int N = Lx*Ly; // 250x40 right now
	int left_Border = 0;
	int right_Border = 0.1*Ly; // Ly*(1.0 / 100.0);  // gives correct integer
	double ratio = (right_Border - left_Border) / (1.0*Ly);
	double Temp_KLS = 3.2;  // kT/J actually
	double Temp_DDS = 5.0;  // kT/J actually
	int timeWindow = 100;
	int totalMCS = 1000; // Total number of Monte Carlo steps per single run
	int runsNumber = 1;

	// ==================================================
	// Memory allocation
	// ==================================================

	int** SpinSys = new int*[Lx];
	int** ZeroTimeSpSys = new int*[Lx];
	// Correlation function for the whole system
	double* Corr = new double[totalMCS];
	double* CorrTot = new double[totalMCS];
	// Correlation function for the particular region
	double* Corr_DDS = new double[totalMCS];
	double* CorrTot_DDS = new double[totalMCS];
	double* Corr_KLS = new double[totalMCS];
	double* CorrTot_KLS = new double[totalMCS];
	for (i = 0; i < Lx; i++) {
		SpinSys[i] = new int[Ly];
		ZeroTimeSpSys[i] = new int[Ly];
	}
	Corr[0] = N / 2 - N*N / 4;
	CorrTot[0] = N / 2 - N*N / 4;

	// ==========================================================================================
	// Precompute all possible probabilities for spins exchangeand some values for optimization
	// ==========================================================================================

	int tableSize = 6, SpinSum = -12;
	double* ExchProb_KLS = new double[tableSize];
	double* ExchProb_DDS = new double[tableSize];
	for (int n = 0; n < 7; n++) {
		ExchProb_KLS[n] = exp(-SpinSum / Temp_KLS);
		ExchProb_DDS[n] = exp(-SpinSum / Temp_DDS);
		if (ExchProb_KLS[n] > 1) {
			ExchProb_KLS[n] = 1;
		}
		if (ExchProb_DDS[n] > 1) {
			ExchProb_DDS[n] = 1;
		}
		SpinSum += 4;
	}
	double size = 1.0*Lx*Ly;

	// ==================================================
	// Random numbers generator creating
	// ==================================================

	random_device rd{};
	mt19937 RNG{ rd() };
	uniform_int_distribution<int> Latt{ 0, Lx*Ly - 1 }; // Here I mapped 2d array to 1d
	uniform_real_distribution<double> Rand{ 0.0, 1.0 };

	// ==================================================
	// Simulation
	// ==================================================

	for (int iwalk = 0; iwalk < runsNumber; iwalk++) {

		//// Clearing all lattice cells from the last run ////
		for (i = 0; i < Lx; i++) {
			for (j = 0; j < Ly; j++) {
				SpinSys[i][j] = 0;
				ZeroTimeSpSys[i][j] = 0;
			}
		}

		//// Clearing array from the last run  ////
		for (dt = 1; dt < totalMCS; dt++) {
			Corr[dt] = 0;
			Corr_DDS[dt] = 0;
			Corr_KLS[dt] = 0;
		}

		//// Filling the lattice with spin particles in random fashion ////
		int positSpins = 0;
		for (n = 0; n < N; n++) {
			// Pick random cell
			dice = Latt(RNG);
			i = dice / Ly; j = dice - i*Ly;
			if (SpinSys[i][j] == 0) {// as soon as positSpins reachees N/2 value - fill the reso of the lattice with negative spins
				if (positSpins < 0.5*N) {
					SpinSys[i][j] = 1;
					positSpins++;
				}else {
					SpinSys[i][j] = -1;
				}
			}else {
				n--;
			}
		}

		//// Single run ////
		for (int istep = 1; istep < totalMCS; istep++) {
			for (int moveAttempt = 0; moveAttempt < N; moveAttempt++) {
				dice = Latt(RNG); // Picking the random spin in the array
				i = dice / Ly; j = dice - i*Ly;
				if (SpinSys[i][j] == 1) { // Work only with positive spins
					iPrev = i == 0 ? Lx - 1 : i - 1;	// if i == 0, then iPrev = Lx-1, else iPrev = i-1;
					if (i == 0) {
						iPPrev = Lx - 2;
					}
					else if (i == 1) {
						iPPrev = Lx - 1;
					}
					else {
						iPPrev = i - 2;
					}
					iNext = i == Lx - 1 ? 0 : i + 1;
					if (i == Lx - 1) {
						iNNext = 1;
					}
					else if (i == Lx - 2) {
						iNNext = 0;
					}
					else {
						iNNext = i + 2;
					}
					jPrev = j == 0 ? Ly - 1 : j - 1;
					jNext = j == Ly - 1 ? 0 : j + 1;

					r = Rand(RNG);
					if (r < 0.5) {	// moving to the right
						if (SpinSys[i][jNext] == -1) { // If atom the right neighbor cell is empty
							SpinSys[i][jNext] = 1;
							SpinSys[i][j] = -1;
						}
					}
					if (0.5 < r&& r < 0.75) {// moving up
						E1 = SpinSys[iPPrev][j] - SpinSys[iNext][j] - SpinSys[i][jPrev] - SpinSys[i][jNext] + SpinSys[iPrev][jPrev] + SpinSys[iPrev][jNext];
						delta = -2 * E1; // E2 = -1 * E1; delta = E2 - E1;
						if (j >= left_Border && j < right_Border) {
							W = ExchProb_DDS[(delta + 12) / 4];
						}
						else {
							W = ExchProb_KLS[(delta + 12) / 4];
						}
						w = Rand(RNG);
						if (w < W) { // Moving spin up by exchange
							SpinSys[i][j] = -1;
							SpinSys[iPrev][j] = 1;
						}
					}
					if (0.75 < r && r < 1.0) {// moving down	
						E1 = -SpinSys[iPrev][j] + SpinSys[iNNext][j] - SpinSys[i][jPrev] - SpinSys[i][jNext] + SpinSys[iNext][jPrev] + SpinSys[iNext][jNext];
						delta = -2 * E1; // E2 = -1 * E1; delta = E2 - E1;
						if (j >= left_Border && j < right_Border) {
							W = ExchProb_DDS[(delta + 12) / 4];
						}
						else {
							W = ExchProb_KLS[(delta + 12) / 4];
						}
						w = Rand(RNG);
						if (w < W) { // Moving spin down by exchange
							SpinSys[i][j] = -1;
							SpinSys[iNext][j] = 1;
						}
					}

				}
			}
			//// Computing correlation function ////
			if (istep == timeWindow) {
				for (i = 0; i < Lx; i++) {
					for (j = 0; j < Ly; j++) {
						ZeroTimeSpSys[i][j] = 0.5*(1 + SpinSys[i][j]); // Mapping  +1/-1 spin system to the 0/1 density system
					}
				}
			}
			// Correlation funciton for the whole region
			if (istep > timeWindow) {
				int Sum = 0;
				int Sum_DDS = 0;
				int Sum_KLS = 0;
				for (i = 0; i < Lx; i++) {
					for (j = 0; j < Ly; j++) {
						if (SpinSys[i][j] == -1) {
							n = 0;
						}
						else {
							n = 1;
						}
						Sum += n * ZeroTimeSpSys[i][j];
						if (j < right_Border) {
							Sum_DDS += n * ZeroTimeSpSys[i][j];
						}
						else {
							Sum_KLS += n * ZeroTimeSpSys[i][j];
						}
					}
				}
				Corr_DDS[istep] = 1.0*Sum_DDS / size - 0.25*ratio;
				Corr_KLS[istep] = 1.0*Sum_KLS / size - 0.25*(1.0 - ratio);
				Corr[istep] = (SpinSum / size)*1.0 - 0.25;
			}
		}

		//// Computing parameters after single run  ////
		for (dt = 1; dt < totalMCS; dt++) {
			CorrTot[dt] += Corr[dt];
			CorrTot_DDS[dt] += Corr_DDS[dt];
			CorrTot_KLS[dt] += Corr_KLS[dt];
		}
	}

	// ==================================================
	// Simulation results output
	// ==================================================

	string outFile = "output.txt";
	ofstream out_stream;
	out_stream.open(outFile.c_str());
	//cout << "MCS" << " " << "t/s" << " " << "S(s,t)" << " " << "S(s,t)*s^0.5" << " " << "S(s,t)*s" << " " << "S_DDS(s,t)" << " " << "S_DDS(s,t)*s^0.5" << " " << "S_DDS(s,t)*s" << " " << "S_KLS(s,t)" << " " << "S_KLS(s,t)*s^0.5" << " " << "S_KLS(s,t)*s" << endl;
	for (dt = timeWindow + 1; dt < totalMCS; dt++) {
		CorrTot[dt] = CorrTot[dt] / runsNumber;
		CorrTot_KLS[dt] = CorrTot_KLS[dt] / runsNumber;
		CorrTot_DDS[dt] = CorrTot_DDS[dt] / runsNumber;
		out_stream << dt << " " << 1.0*dt / (1.0*timeWindow) << "	" << CorrTot[dt] << " " << CorrTot[dt] * (pow(timeWindow, 0.5)) << " " << CorrTot[dt] * (pow(timeWindow, 1.0)) << " " << CorrTot_DDS[dt] << " " << CorrTot_DDS[dt] * pow(timeWindow, 0.5) << " " << CorrTot_DDS[dt] * pow(timeWindow, 1.0) << " " << CorrTot_KLS[dt] << " " << CorrTot_KLS[dt] * pow(timeWindow, 0.5) << " " << CorrTot_KLS[dt] * pow(timeWindow, 1.0) << endl;
	}
	out_stream.close();

	// ==================================================
	// Memory deallocation
	// ==================================================

	delete[] SpinSys;
	delete[] ZeroTimeSpSys;
	delete[] Corr;
	delete[] CorrTot;
	delete[] Corr_DDS;
	delete[] Corr_KLS;
	delete[] CorrTot_DDS;
	delete[] CorrTot_KLS;

	return 0;

}
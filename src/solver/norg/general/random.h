#pragma once


#include <random>
#include <iostream>
#include <iomanip>
#ifdef _MSC_VER
#include "..\\randomc\\randomc.h"
#include "..\\randomc\\sfmt.h"
#else
#include "../randomc/randomc.h"
#include "../randomc/sfmt.h"
#endif
// uniform unit random number [0, 1)
class UURandom {
private:
	virtual double uniform_unit_random() = 0;
	unsigned parallel_seed(int myid, unsigned seed) {
		// you can change global_shift to get a new set of parallel seeds
		CRandomSFMT rng(seed);
		for (int i = 0; i < myid; i++) {
			rng.IRandom(0, 0x7fffffff);
		}
		return rng.IRandom(0, 0x7fffffff);
	}
public:
	virtual ~UURandom() {}
	virtual void init_random(int seed) = 0;
	void init_random_parallel(int myid, unsigned seed) {
		unsigned seed_parallel = parallel_seed(myid, seed);
		// std::cout << "random seed = " << std::setw(10) << seed_parallel << " with myid = " << myid << std::endl;
		init_random(seed_parallel);
	}


	double operator()() {
		return uniform_unit_random();
	}
	void operator()(VecReal& v) {
		for_Int(i, 0, v.size()) {
			v[i] = (*this)();
		}
	}
	void operator()(VecCmplx& v) {
		for_Int(i, 0, v.size()) {
			v[i] = Cmplx((*this)(), (*this)());
		}
	}
};



// uniform unit random number [0, 1), Mersenne
class UURandomMersenne : public UURandom {
private:
	CRandomMersenne rng;
	double uniform_unit_random() { return rng.Random(); }
public:
	UURandomMersenne(int myid = 0, int seed = 1) : rng(seed) { init_random_parallel(myid, unsigned(seed)); }
	void init_random(int seed) { rng.RandomInit(seed); }
};


// uniform unit random number [0, 1), SFMT
class UURandomSFMT : public UURandom {
private:
	CRandomSFMT rng;
	double uniform_unit_random() { return rng.Random(); }
public:
	UURandomSFMT(int myid = 0, int seed = 1) : rng(seed) { init_random_parallel(myid, unsigned(seed)); }
	void init_random(int seed) { rng.RandomInit(seed); }
};





#ifdef _MSC_VER
// uniform unit random number [0, 1), STL
class UURandomSTL : public UURandom {
private:
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution;
public:
	UURandomSTL(int myid = 0, int seed = 1) : generator(seed), distribution(0., 1.) { init_random_parallel(myid, unsigned(seed)); }
	void init_random(int seed) { generator.seed(seed); }
	double uniform_unit_random() { return distribution(generator); }
};
#endif
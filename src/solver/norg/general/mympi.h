/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China) date 2013 - 2017
	Jia-Ming Wang (jmw@ruc.edu.cn, RUC, China) date 2020 - 2022
*/

#ifndef _MYMPI_H_
#define _MYMPI_H_

#include "stdpfx.h"
#include "mpi.h"
#include "vec.h"
#include "myomp.h"
#include "safebool.h"

class VecPartition {
private:
	Int np;		// number of processes
	Int id;		// myid for this process
	Int n;		// [0, 1, ..., n - 1] is what you want to partition
	VecInt _bgn;
	VecInt _end;
	VecInt _len;
	VecInt _dsp;
public:
	VecPartition(Int np_i, Int id_i, Int n_i): np(np_i), id(id_i), n(n_i), _bgn(np), _end(np), _len(np) {
		if (np < 1) ERR(NAV2(np, n));
		Int k = n / np;
		Int j = n % np;
		Int m = 0;
		for_Int (i, 0, np) {
			_bgn[i] = m;
			m += k + (i < j);
			_len[i] = m - _bgn[i];
			_end[i] = m;
		}
		if (m != n) ERR(NAV3(np, n, m));
		_dsp = _bgn;
	}
	Int total_size() const { return n; }
	Int bgn() const { return _bgn[id]; }
	Int end() const { return _end[id]; }
	Int len() const { return _len[id]; }
	Int dsp() const { return _dsp[id]; }
	Int bgn(const Int &id_i) const { return _bgn[id_i]; }
	Int end(const Int &id_i) const { return _end[id_i]; }
	Int len(const Int &id_i) const { return _len[id_i]; }
	Int dsp(const Int &id_i) const { return _dsp[id_i]; }
	const Int* bgn_p() const { return _bgn.p(); }
	const Int* end_p() const { return _end.p(); }
	const Int* len_p() const { return _len.p(); }
	const Int* dsp_p() const { return _dsp.p(); }
};

class MyMpi: public safe_bool<MyMpi> {
private:
	typedef std::complex<double> _Cmplx;
	MPI_Comm comm;
	int nprc;
	int myid;
	int mstr;
private:
	int Send(char* buf, int count, int dest, int tag) const {
		return MPI_Send(buf, count, MPI_CHAR, dest, tag, comm);
	}
	int Send(int* buf, int count, int dest, int tag) const {
		return MPI_Send(buf, count, MPI_INT, dest, tag, comm);
	}
	int Send(double* buf, int count, int dest, int tag) const {
		return MPI_Send(buf, count, MPI_DOUBLE, dest, tag, comm);
	}
	int Send(_Cmplx* buf, int count, int dest, int tag) const {
		return MPI_Send(buf, count, MPI_DOUBLE_COMPLEX, dest, tag, comm);
	}
	int Recv(char* buf, int count, int source, int tag, MPI_Status* status) const {
		return MPI_Recv(buf, count, MPI_CHAR, source, tag, comm, status);
	}
	int Recv(int* buf, int count, int source, int tag, MPI_Status* status) const {
		return MPI_Recv(buf, count, MPI_INT, source, tag, comm, status);
	}
	int Recv(double* buf, int count, int source, int tag, MPI_Status* status) const {
		return MPI_Recv(buf, count, MPI_DOUBLE, source, tag, comm, status);
	}
	int Recv(_Cmplx* buf, int count, int source, int tag, MPI_Status* status) const {
		return MPI_Recv(buf, count, MPI_DOUBLE_COMPLEX, source, tag, comm, status);
	}
	int Bcast(char *buff, int count) const {
		return MPI_Bcast(buff, count, MPI_CHAR, mstr, comm);
	}
	int Bcast(int *buff, int count) const {
		return MPI_Bcast(buff, count, MPI_INT, mstr, comm);
	}
	int Bcast(double *buff, int count) const {
		return MPI_Bcast(buff, count, MPI_DOUBLE, mstr, comm);
	}
	int Bcast(_Cmplx *buff, int count) const {
		return MPI_Bcast(buff, count, MPI_DOUBLE_COMPLEX, mstr, comm);
	}

	int Reduce(const char *sendbuf, char *recvbuf, int sendcount) const {
		return MPI_Reduce(sendbuf, recvbuf, sendcount, MPI_CHAR, MPI_SUM, mstr, comm);
	}
	int Reduce(const int *sendbuf, int *recvbuf, int sendcount) const {
		return MPI_Reduce(sendbuf, recvbuf, sendcount, MPI_INT, MPI_SUM, mstr, comm);
	}
	int Reduce(const double *sendbuf, double *recvbuf, int sendcount) const {
		return MPI_Reduce(sendbuf, recvbuf, sendcount, MPI_DOUBLE, MPI_SUM, mstr, comm);
	}
	int Reduce(const _Cmplx *sendbuf, _Cmplx *recvbuf, int sendcount) const {
		return MPI_Reduce(sendbuf, recvbuf, sendcount, MPI_DOUBLE_COMPLEX, MPI_SUM, mstr, comm);
	}
	int Allreduce(const char *sendbuf, char *recvbuf, int sendcount) const {
		return MPI_Allreduce(sendbuf, recvbuf, sendcount, MPI_CHAR, MPI_SUM, comm);
	}
	int Allreduce(const int *sendbuf, int *recvbuf, int sendcount) const {
		return MPI_Allreduce(sendbuf, recvbuf, sendcount, MPI_INT, MPI_SUM, comm);
	}
	int Allreduce(const double *sendbuf, double *recvbuf, int sendcount) const {
		return MPI_Allreduce(sendbuf, recvbuf, sendcount, MPI_DOUBLE, MPI_SUM, comm);
	}
	int Allreduce(const _Cmplx *sendbuf, _Cmplx *recvbuf, int sendcount) const {
		return MPI_Allreduce(sendbuf, recvbuf, sendcount, MPI_DOUBLE_COMPLEX, MPI_SUM, comm);
	}
	int Gatherv(const double *sendbuf, int sendcount, double *recvbuf, const int *recvcount, const int *displs) const {
		return MPI_Gatherv(sendbuf, sendcount, MPI_DOUBLE, recvbuf, recvcount, displs, MPI_DOUBLE, mstr, comm);
	}
	int Gatherv(const _Cmplx *sendbuf, int sendcount, _Cmplx *recvbuf, const int *recvcount, const int *displs) const {
		return MPI_Gatherv(sendbuf, sendcount, MPI_DOUBLE_COMPLEX, recvbuf, recvcount, displs, MPI_DOUBLE_COMPLEX, mstr, comm);
	}
	int Gather(const double *sendbuf, int sendcount, double *recvbuf, int recvcount) const {
		return MPI_Gather(sendbuf, sendcount, MPI_DOUBLE, recvbuf, recvcount, MPI_DOUBLE, mstr, comm);
	}
	int Gather(const _Cmplx *sendbuf, int sendcount, _Cmplx *recvbuf, int recvcount) const {
		return MPI_Gather(sendbuf, sendcount, MPI_DOUBLE_COMPLEX, recvbuf, recvcount, MPI_DOUBLE_COMPLEX, mstr, comm);
	}
	int Allgatherv(const int* sendbuf, int sendcount, int* recvbuf, const int* recvcount, const int* displs) const {
		return MPI_Allgatherv(sendbuf, sendcount, MPI_INT, recvbuf, recvcount, displs, MPI_INT, comm);
	}
	int Allgatherv(const double* sendbuf, int sendcount, double* recvbuf, const int* recvcount, const int* displs) const {
		return MPI_Allgatherv(sendbuf, sendcount, MPI_DOUBLE, recvbuf, recvcount, displs, MPI_DOUBLE, comm);
	}
	int Allgatherv(const _Cmplx* sendbuf, int sendcount, _Cmplx* recvbuf, const int* recvcount, const int* displs) const {
		return MPI_Allgatherv(sendbuf, sendcount, MPI_DOUBLE_COMPLEX, recvbuf, recvcount, displs, MPI_DOUBLE_COMPLEX, comm);
	}
	int Allgather(const char *sendbuf, int sendcount, char *recvbuf, int recvcount) const {
		return MPI_Allgather(sendbuf, sendcount, MPI_CHAR, recvbuf, recvcount, MPI_CHAR, comm);
	}
	int Allgather(const double *sendbuf, int sendcount, double *recvbuf, int recvcount) const {
		return MPI_Allgather(sendbuf, sendcount, MPI_DOUBLE, recvbuf, recvcount, MPI_DOUBLE, comm);
	}
	int Allgather(const _Cmplx* sendbuf, int sendcount, _Cmplx* recvbuf, int recvcount) const {
		return MPI_Allgather(sendbuf, sendcount, MPI_DOUBLE, recvbuf, recvcount, MPI_DOUBLE, comm);
	}
public:
	MyMpi(MPI_Comm comm_i = MPI_COMM_WORLD, int mstr_i = 0) : comm(comm_i), mstr(mstr_i) {
		MPI_Comm_size(comm, &nprc);
		MPI_Comm_rank(comm, &myid);
		if (myid < 0 || myid >= nprc || mstr < 0 || mstr >= nprc) ERR(NAV4(comm, nprc, myid, mstr));
		display();
	}
	MPI_Comm cm() const { return comm; }
	int np() const { return nprc; }
	int id() const { return myid; }
	int ms() const { return mstr; }
	bool boolean_test() const { return myid == mstr; }
	int barrier() const { return MPI_Barrier(comm); }
	void print_process_names(std::ostream &os = std::cout) const;
	void display() const {
		if (*this) std::cout << "MPI: number of processes = " << np() << std::endl;
		if (*this) std::cout << "MPI: master process id = " << ms() << std::endl;
		print_process_names();
	}
public:
	template <typename T>
	int Send(const T& v, int dest, int tag) const {
		return Send(&v, 1, dest, tag);
	}
	template <typename T>
	int Send(const Vec<T>& v, int dest, int tag) const {
		return Send(v.p(), (int)v.size(), dest, tag);
	}
	template <typename T>
	int Send(const Mat<T>& a, int dest, int tag) const {
		return Send(a.p(), (int)a.size(), dest, tag);
	}
	template <typename T>
	int Recv(T& v, MPI_Status& status, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG) const {
		return Recv(&v, 1, source, tag, &status);
	}
	template <typename T>
	int Recv(Vec<T>& v, MPI_Status& status, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG) const {
		return Recv(v.p(), (int)v.size(), source, tag, &status);
	}
	template <typename T>
	int Recv(Mat<T>& a, MPI_Status& status, int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG) const {
		return Recv(a.p(), (int)a.size(), source, tag, &status);
	}
	template <typename T>
	int Bcast(T &v) const {
		return Bcast(&v, 1);
	}
	template <typename T>
	int Bcast(Vec<T> &v) const {
		return Bcast(v.p(), (int)v.size());
	}
	template <typename T>
	int Bcast(Mat<T>& a) const {
		return Bcast(a.p(), (int)a.size());
	}
	template <typename T>
	T Reduce(const T &sendbuf) const {
		T recvbuf;
		Reduce(&sendbuf, &recvbuf, 1);
		return recvbuf;
	}
	template <typename T>
	Vec<T> Reduce(const Vec<T> &sendbuf) const {
		Vec<T> recvbuf(*this ? sendbuf.size() : 0);
		Reduce(sendbuf.p(), recvbuf.p(), (int)sendbuf.size());
		return *this ? recvbuf : sendbuf;
	}
	template <typename T>
	Mat<T> Reduce(const Mat<T> &sendbuf) const {
		Mat<T> recvbuf(*this ? sendbuf.nrows() : 0, *this ? sendbuf.ncols() : 0);
		Reduce(sendbuf.p(), recvbuf.p(), (int)sendbuf.size());
		return *this ? recvbuf : sendbuf;
	}
	template <typename T>
	T Allreduce(const T &sendbuf) const {
		T recvbuf;
		Allreduce(&sendbuf, &recvbuf, 1);
		return recvbuf;
	}
	template <typename T>
	Vec<T> Allreduce(const Vec<T> &sendbuf) const {
		Vec<T> recvbuf(sendbuf.size());
		Allreduce(sendbuf.p(), recvbuf.p(), (int)sendbuf.size());
		return recvbuf;
	}
	template <typename T>
	Mat<T> Allreduce(const Mat<T>& sendbuf) const {
		Mat<T> recvbuf(sendbuf.nrows(), sendbuf.ncols());
		Allreduce(sendbuf.p(), recvbuf.p(), (int)sendbuf.size());
		return recvbuf;
	}
	template <typename T>
	Vec<T> Gatherv(const Vec<T> &sendbuf, const VecPartition &vp) const {
		Vec<T> recvbuf(*this ? vp.total_size() : 0);
		Gatherv(sendbuf.p(), (int)sendbuf.size(), recvbuf.p(), (const int*)vp.len_p(), (const int*)vp.dsp_p());
		return recvbuf;
	}
	template <typename T>
	Vec<Mat<T>> Gatherv(const Vec<Mat<T>>& vm, Idx mat_nrows, Idx mat_ncols, const VecPartition& vp) const;

	template <typename T>
	Vec<T> Allgatherv(const Vec<T>& sendbuf, const VecPartition& vp) const {
		Vec<T> recvbuf(vp.total_size());
		Allgatherv(sendbuf.p(), (int)sendbuf.size(), recvbuf.p(), (const int*)vp.len_p(), (const int*)vp.dsp_p());
		return recvbuf;
	}
	template <typename T>
	Vec<T> Gather(const T &sendbuf) const {
		Vec<T> recvbuf(*this ? np() : 0);
		Gather(&sendbuf, 1, recvbuf.p(), 1);
		return recvbuf;
	}
	template <typename T>
	Vec<T> Gather(const Vec<T> &sendbuf) const {
		Vec<T> recvbuf((*this ? np() : 0) * sendbuf.size());
		Gather(sendbuf.p(), (int)sendbuf.size(), recvbuf.p(), (int)sendbuf.size());
		return recvbuf;
	}
	template <typename T>
	Vec<T> Allgather(const T &sendbuf) const {
		Vec<T> recvbuf(np());
		Allgather(&sendbuf, 1, recvbuf.p(), 1);
		return recvbuf;
	}
	template <typename T>
	Vec<T> Allgather(const Vec<T> &sendbuf) const {
		Vec<T> recvbuf(np() * sendbuf.size());
		Allgather(sendbuf.p(), (int)sendbuf.size(), recvbuf.p(), (int)sendbuf.size());
		return recvbuf;
	}
	int Allreduce_MIN(const double* sendbuf, double* recvbuf, int sendcount) const {
		return MPI_Allreduce(sendbuf, recvbuf, sendcount, MPI_DOUBLE, MPI_MIN, comm);
	}
};

template <typename T>
Vec<Mat<T>> MyMpi::Gatherv(const Vec<Mat<T>>& vm, Idx mat_nrows, Idx mat_ncols, const VecPartition& vp) const
{
	Mat<Vec<T>> sendbuf(mat_nrows, mat_ncols, Vec<T>(vm.size()));
	for_Idx(i, 0, mat_nrows) {
		for_Idx(j, 0, mat_ncols) {
			for_Idx(k, 0, vm.size()) {
				sendbuf[i][j][k] = vm[k][i][j];
			}
		}
	}

	Mat<Vec<T>> recvbuf(mat_nrows, mat_ncols);
	for_Idx(i, 0, mat_nrows) {
		for_Idx(j, 0, mat_ncols) {
			recvbuf[i][j] = Gatherv(sendbuf[i][j], vp);
		}
	}

	Idx len = (*this ? vp.total_size() : 0);
	Vec<Mat<T>> res(len, Mat<T>(mat_nrows, mat_ncols));
	for_Idx(i, 0, mat_nrows) {
		for_Idx(j, 0, mat_ncols) {
			for_Idx(k, 0, len) {
				res[k][i][j] = recvbuf[i][j][k];
			}
		}
	}
	return res;
}

inline void use_omp(const MyMpi &mm, Int nthreads = 1)
{
#ifdef _OPENMP
	omp_set_num_threads(nthreads);
	int omp_num_threads;
#pragma omp parallel
	omp_num_threads = omp_get_num_threads();
	if (mm) std::cout << "OpenMP: number of procs = " << omp_get_num_procs() << std::endl;
	if (mm) std::cout << "OpenMP: omp_num_threads = " << omp_num_threads << std::endl;
#endif
}

inline void use_mkl(const MyMpi &mm, Int nthreads = 1)
{
	//mkl_set_num_threads(nthreads);
	//mkl_set_dynamic;
	if (mm) std::cout << "MKL: mkl_num_threads = " << mkl_get_max_threads() << std::endl;
}

#endif /* _MYMPI_H_ */

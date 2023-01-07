/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2013 - 2017
*/

#include "mympi.h"

void MyMpi::print_process_names(std::ostream &os) const
{
	char tmp[MPI_MAX_PROCESSOR_NAME];
	char pnp[MPI_MAX_PROCESSOR_NAME];
	int pnl;
	MPI_Status stts;

	MPI_Get_processor_name(pnp, &pnl);
	if (*this) {
		int len_before_pid = TOSTRLEN("Host Name");
		if (len_before_pid < TOSTRLEN(pnp)) len_before_pid = TOSTRLEN(pnp);
		for_Int (pid, 0, this->np()) if (pid != this->id()) {
			MPI_Recv(tmp, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, pid, 0, this->cm(), &stts);
			if (len_before_pid < TOSTRLEN(tmp)) len_before_pid = TOSTRLEN(tmp);
		}
		len_before_pid += 4;

		int w_myid = TOSTRLEN(this->np() - 1) + 1;
		if (w_myid < TOSTRLEN("PID")) w_myid = TOSTRLEN("PID");

		os << left_justify("Host Name", len_before_pid);
		os << rght_justify("PID", w_myid) << std::endl;

		Str pnps(pnp);
		os << left_justify(pnps, len_before_pid);
		os << rght_justify(this->id(), w_myid);
		for_Int (pid, 0, this->np()) if (pid != this->id()) {
			MPI_Recv(tmp, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, pid, 1, this->cm(), &stts);
			Str tmps(tmp);
			if (pnps != tmps) {
				pnps = tmps;
				os << std::endl;
				os << left_justify(pnps, len_before_pid);
				os << rght_justify(pid, w_myid);;
			} else {
				os << rght_justify(pid, w_myid);
			}
		}
		os << std::endl;
	} else {
		MPI_Send(pnp, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, this->ms(), 0, this->cm());
		MPI_Send(pnp, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, this->ms(), 1, this->cm());
	}
}

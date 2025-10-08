#include <stdio.h>
#include "zlib.h"
#include "mpi.h"

typedef struct {
  uLong checksum;
  z_off_t len;
} checktype;
  
void dCRC_red(checktype *invec, checktype *inoutvec, int *len,
              MPI_Datatype *datatype) {
  int i;
  checktype combi;

  for (i=0; i< *len; ++i) {
    combi.len = invec->len + inoutvec->len;
    combi.checksum = crc32_combine(invec->checksum,
                                   inoutvec->checksum,
                                   inoutvec->len);
    *inoutvec = combi;
    invec++; inoutvec++;
  }
}

void distCRC(const void *dat, long length, MPI_Fint comm, char buf[]) {
  MPI_Comm world;
  MPI_Op myOp;
  MPI_Datatype cstype;

  checktype cs, res;

  world = MPI_Comm_f2c(comm);
  MPI_Type_contiguous(2, MPI_LONG, &cstype);
  MPI_Type_commit(&cstype);
  MPI_Op_create(dCRC_red, 0, &myOp);

  cs.len = (z_off_t) length;
  cs.checksum = crc32(0L, dat, length);

  MPI_Allreduce(&cs, &res, 1, cstype, myOp, world);

  MPI_Op_free(&myOp);
  MPI_Type_free(&cstype);

  sprintf(buf, "%8.8x", res.checksum);
}

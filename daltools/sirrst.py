#!/usr/bin/env python
import numpy
import sys
import os
import tarfile
import tempfile
import shutil
from util import blocked, full
from fortran_binary import FortranBinary
from .basinfo import BasInfo

class SiriusRestart(object):
    def __init__(self, name="SIRIUS.RST", tgz=None):
        self.name = name
        if tgz is not None:
            tmp = tempfile.mkdtemp()
            tarfile.open(tgz, 'r:gz').extractall(
                path=tmp
                )
            self.name = os.path.join(tmp, "SIRIUS.RST")
        self.basinfo = BasInfo(self.name)
        #self.close()
        #if tgz:
        #    shutil.rmtree(tmp)
        self._cmo = None
        self._ci = None

    def __str__(self):
        retstr=""
        retstr+="\nNSYM : " + str(self.basinfo.nsym)
        retstr+="\nNBAS : " + str(self.basinfo.nbas)
        retstr+="\nNBAST: " + str(self.basinfo.nbast)
        retstr+="\nNORB : " + str(self.basinfo.norb)
        retstr+="\nNORBT: " + str(self.basinfo.norbt)
        retstr+="\nNRHF : " + str(self.basinfo.nrhf)
        retstr+="\nIOPRHF:" + str(self.basinfo.ioprhf)
        retstr+="\nCMO:   " + str(self.cmo)
        return retstr

    def get_rhf_density(self):
        densities = blocked.BlockDiagonalMatrix(self.basinfo.nbas, self.basinfo.nbas)
        for dens, occ, cmo in zip(densities, self.basinfo.nrhf, self.cmo):
            dens += 2*cmo[:, :occ]*cmo[:, :occ].T
        return densities.unblock()

    def get_occ_density(self, occnum):
        densities = blocked.BlockDiagonalMatrix(self.basinfo.nbas,
                self.basinfo.nbas)
        for dens, occ, cmo in zip(densities, occnum, self.cmo):
            for i, ni in enumerate(occ):
                dens += ni*numpy.outer(cmo[:, i], cmo[:, i])
        return densities.unblock()

    @property
    def cmo(self):
        if self._cmo is None:
            fb = FortranBinary(self.name)
            fb.find("NEWORB")
            cmo_rec=fb.next()
            assert cmo_rec.reclen//8 == numpy.dot(self.basinfo.nbas, self.basinfo.norb)
            ncmot4 = max(self.basinfo.ncmot,4)
            n=0
            cmo=blocked.BlockDiagonalMatrix(self.basinfo.nbas, self.basinfo.norb)
            for isym in range(self.basinfo.nsym):
                cmoi = numpy.array(cmo_rec.read(self.basinfo.nbas[isym]*self.basinfo.norb[isym],'d')
                       ).reshape((self.basinfo.nbas[isym], self.basinfo.norb[isym]), order='F')
                cmo.subblock[isym] = cmoi.view(full.matrix)
            self._cmo = cmo
            fb.close()
        return self._cmo
            
    @property
    def ci(self):
        if self._ci is None:
            fb = FortranBinary(self.name)
            fb.find('STARTVEC')
            ci_record = fb.next()
            fb.close()
            self._ci = numpy.array(ci_record.read(ci_record.reclen//8, 'd'))       
        return self._ci
            


def main():

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('SIRIUS_RST')

    args = parser.parse_args()

    rst=SiriusRestart(args.SIRIUS_RST)
    print(rst)

if __name__ == "__main__":#pragma: no cover
    sys.exit(main())

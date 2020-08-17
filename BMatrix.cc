#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

#include "K_matrix_base.h"
#include "box_quant.h"

namespace py = pybind11;


  /* A valid K matrix parametrization requires two functions which behave
   * according to
   *
   *  double calculate(uint Jtimestwo,
   *                   uint Lp, uint Sptimestwo, uint chanp,
   *                   uint L, uint Stimestwo, uint chan,
   *                   double Ecm_over_mref) const;
   *
   *  bool isZero(uint Jtimestwo,
   *              uint Lp, uint Sptimestwo, uint chanp,
   *              uint L, uint Stimestwo, uint chan) const;
   *
   *
   * This class is a bit stupid because it only forwards the function calls
   * to the two std::functions passed in.
   */

class KMatrix : public KtildeMatrixBase {

  public:

    typedef std::vector<std::function<double(double)>> qSqFctsType;
    typedef std::function<double(uint,uint,uint,uint,uint,uint,uint,double,const qSqFctsType&)> calcFuncType;
    typedef std::function<bool(uint,uint,uint,uint,uint,uint,uint)> zeroFuncType;

    KMatrix(const calcFuncType& _calculate, zeroFuncType _isZero) : calcFunc(_calculate), isZeroFunc(_isZero) {};

    double calculate(uint Jtimestwo,
                     uint Lp, uint Sptimestwo, uint chanp,
                     uint L, uint Stimestwo, uint chan,
                     double Ecm_over_mref,
		     const qSqFctsType& qSqFcts) const { return calcFunc(Jtimestwo, Lp, Sptimestwo, chanp, L, Stimestwo, chan, Ecm_over_mref, qSqFcts); }

    bool isZero(uint Jtimestwo,
                uint Lp, uint Sptimestwo, uint chanp,
                uint L, uint Stimestwo, uint chan) const { return isZeroFunc(Jtimestwo, Lp, Sptimestwo, chanp, L, Stimestwo, chan); }

    uint getNumberOfParameters() const { return 1; }


  private:

    calcFuncType calcFunc;
    zeroFuncType isZeroFunc;

};


class FastBoxQuantization {

  BoxQuantization BQ;
  std::vector<double> elabList, ecmList;
  std::vector<ComplexHermitianMatrix> bMatList;

  public:

    FastBoxQuantization(const std::string& mom_ray, uint mom_int_sq,
		        const std::string& lgirrep,
		        const std::vector<DecayChannelInfo>& chan_infos,
		        const std::vector<uint> lmaxes,
		        KMatrix& K,
		        bool KInvMode)
      	: BQ(mom_ray, mom_int_sq, lgirrep, chan_infos, lmaxes, &K, KInvMode)
		{}

    void initialize(std::vector<double>& _elabList,
	       std::vector<double>& _refMassLList,
	       std::list<std::pair<std::vector<double>,std::vector<double>>>& _shMassList
	      ) {

      elabList = _elabList;
      ecmList.reserve(elabList.size());

      for (unsigned int iS = 0; iS < _elabList.size(); iS++) {

	BQ.setRefMassL(_refMassLList[iS]);

	unsigned int iChan = 0;
	for (auto chanIt = _shMassList.begin(); chanIt != _shMassList.end(); ++chanIt, ++iChan)
	  BQ.setMassesOverRef(iChan, chanIt->first[iS], chanIt->second[iS]);

	ecmList.push_back(BQ.getEcmOverMrefFromElab(_elabList[iS]));

	ComplexHermitianMatrix bMat;
	BQ.getBoxMatrixFromEcm(ecmList[iS], bMat);
	bMatList.push_back(bMat);

      } // loop over samples
    }

    double getOmegaForSample(double mu, unsigned int iS) {
      	// get Kinv or K for Ecm

      RealSymmetricMatrix kMat;
      BQ.getKtildeOrInverseFromEcm(ecmList[iS], kMat);

      return BQ.getOmega(mu, kMat, bMatList[iS]);
    }


    std::vector<double> getOmegaList(double mu) {
      std::vector<double> ret;

      unsigned int nSmpls = ecmList.size();
      ret.reserve(nSmpls);

      for (unsigned int iS = 0; iS < nSmpls; ++iS)
	ret.push_back(getOmegaForSample(mu, iS));

      return ret;
    }

    std::vector<std::complex<double>> getBoxMatrixElementList(unsigned int iRow = 0, unsigned int iCol = 0) {
      std::vector<std::complex<double>> ret;

      unsigned int nSmpls = bMatList.size();
      ret.reserve(nSmpls);

      for (unsigned int iS = 0; iS < nSmpls; ++iS)
	ret.push_back(bMatList[iS](iRow, iCol));

      return ret;
    }

    std::vector<double> getEcmOverMrefList() const { return ecmList; }
    std::vector<double> getElabOverMrefList() const { return elabList; }
};


PYBIND11_MODULE(BMat, m) {
  py::class_<KMatrix>(m, "KMatrix")
    .def(py::init<const KMatrix::calcFuncType&, const KMatrix::zeroFuncType&>())
    .def("calculate", &KMatrix::calculate)
    .def("getName", &KMatrix::isZero);

  py::class_<DecayChannelInfo>(m, "DecayChannelInfo")
    .def(py::init<const std::string&, const std::string&, uint,  uint, bool, bool>());


  	//	BoxQuantization binding of the original class
  py::class_<BoxQuantization>(m, "BoxQuantization")
    .def(py::init([](const std::string& mom_ray, uint mom_int_sq,
                     const std::string& lgirrep,
                     const std::vector<DecayChannelInfo>& chan_infos,
                     const std::vector<uint> lmaxes,
                     KMatrix& K,
		     bool KInvMode)
		{
			    return new BoxQuantization(mom_ray, mom_int_sq, lgirrep, chan_infos, lmaxes, &K, KInvMode);
		}))
    .def("setRefMassL", &BoxQuantization::setRefMassL)
    .def("setMassesOverRef", &BoxQuantization::setMassesOverRef)
    .def("getBoxMatrixFromEcm", [](BoxQuantization& self, double Ecm_over_mref)
		    {
		      ComplexHermitianMatrix res;
		      self.getBoxMatrixFromEcm(Ecm_over_mref, res);
		      return res(0,0);
		    })
    .def("getBoxMatrixFromElab", [](BoxQuantization& self, double Elab_over_mref)
		    {
		      ComplexHermitianMatrix res;
		      self.getBoxMatrixFromElab(Elab_over_mref, res);
		      return res(0,0);
		    })
    .def("getOmegaFromEcm", &BoxQuantization::getOmegaFromEcm)
    .def("getQcmSqFunctions", &BoxQuantization::getQcmSqFunctions)
    .def("getBoxMatrixFromEcmAndSHMasses", [](BoxQuantization& self, double muVal, double Ecm_over_mref, double refMass, const std::vector<std::pair<double, double>>& shMasses) {
		    self.setRefMassL(refMass);

		    uint iChan = 0;
		    for (auto chanIt = shMasses.begin(); chanIt != shMasses.end(); ++chanIt, ++iChan)
		      self.setMassesOverRef(iChan, chanIt->first, chanIt->second);

		    ComplexHermitianMatrix res;
		    self.getBoxMatrixFromEcm(Ecm_over_mref, res);
		    return res(0,0);
		  });

    	//	FastBoxQuantization with a different workflow
  py::class_<FastBoxQuantization>(m, "FastBoxQuantization")
    .def(py::init([](const std::string& mom_ray, uint mom_int_sq,
                     const std::string& lgirrep,
                     const std::vector<DecayChannelInfo>& chan_infos,
                     const std::vector<uint> lmaxes,
                     KMatrix& K,
		     bool KInvMode)
		{
			    return new FastBoxQuantization(mom_ray, mom_int_sq, lgirrep, chan_infos, lmaxes, K, KInvMode);
		}))
    .def("getOmegaForSample", &FastBoxQuantization::getOmegaForSample)
    .def("getOmegaList", &FastBoxQuantization::getOmegaList)
    .def("getBoxMatrixElementList", &FastBoxQuantization::getBoxMatrixElementList)
    .def("getEcmOverMrefList", &FastBoxQuantization::getEcmOverMrefList)
    .def("getElabOverMrefList", &FastBoxQuantization::getElabOverMrefList)
    .def("initialize", &FastBoxQuantization::initialize);
}

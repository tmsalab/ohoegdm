// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// gen_bijectionvector
arma::vec gen_bijectionvector(unsigned int K, unsigned int M);
RcppExport SEXP _ohoegdm_gen_bijectionvector(SEXP KSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(gen_bijectionvector(K, M));
    return rcpp_result_gen;
END_RCPP
}
// inv_gen_bijectionvector
arma::vec inv_gen_bijectionvector(unsigned int K, unsigned int M, double CL);
RcppExport SEXP _ohoegdm_inv_gen_bijectionvector(SEXP KSEXP, SEXP MSEXP, SEXP CLSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< double >::type CL(CLSEXP);
    rcpp_result_gen = Rcpp::wrap(inv_gen_bijectionvector(K, M, CL));
    return rcpp_result_gen;
END_RCPP
}
// CL_gen_invbijection_table
arma::mat CL_gen_invbijection_table(unsigned int K, unsigned int M, unsigned int nClass);
RcppExport SEXP _ohoegdm_CL_gen_invbijection_table(SEXP KSEXP, SEXP MSEXP, SEXP nClassSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    rcpp_result_gen = Rcpp::wrap(CL_gen_invbijection_table(K, M, nClass));
    return rcpp_result_gen;
END_RCPP
}
// ETAmat
arma::mat ETAmat(unsigned int K, unsigned int J, unsigned int M, const arma::mat& Q);
RcppExport SEXP _ohoegdm_ETAmat(SEXP KSEXP, SEXP JSEXP, SEXP MSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(ETAmat(K, J, M, Q));
    return rcpp_result_gen;
END_RCPP
}
// rTruncNorm_lb
double rTruncNorm_lb(double mean, double sd, double b_lb);
RcppExport SEXP _ohoegdm_rTruncNorm_lb(SEXP meanSEXP, SEXP sdSEXP, SEXP b_lbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< double >::type b_lb(b_lbSEXP);
    rcpp_result_gen = Rcpp::wrap(rTruncNorm_lb(mean, sd, b_lb));
    return rcpp_result_gen;
END_RCPP
}
// rTruncNorm
double rTruncNorm(double mean, double sd, double w, const arma::vec& ps);
RcppExport SEXP _ohoegdm_rTruncNorm(SEXP meanSEXP, SEXP sdSEXP, SEXP wSEXP, SEXP psSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ps(psSEXP);
    rcpp_result_gen = Rcpp::wrap(rTruncNorm(mean, sd, w, ps));
    return rcpp_result_gen;
END_RCPP
}
// random_Q
arma::mat random_Q(unsigned int J, unsigned int K);
RcppExport SEXP _ohoegdm_random_Q(SEXP JSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(random_Q(J, K));
    return rcpp_result_gen;
END_RCPP
}
// simSLCM
arma::mat simSLCM(unsigned int N, unsigned int J, unsigned int M, unsigned int nClass, const arma::vec& CLASS, const arma::mat& Atable, const arma::mat& BETA, const arma::mat& TAU);
RcppExport SEXP _ohoegdm_simSLCM(SEXP NSEXP, SEXP JSEXP, SEXP MSEXP, SEXP nClassSEXP, SEXP CLASSSEXP, SEXP AtableSEXP, SEXP BETASEXP, SEXP TAUSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type CLASS(CLASSSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Atable(AtableSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type BETA(BETASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type TAU(TAUSEXP);
    rcpp_result_gen = Rcpp::wrap(simSLCM(N, J, M, nClass, CLASS, Atable, BETA, TAU));
    return rcpp_result_gen;
END_RCPP
}
// BetatoTheta
arma::mat BetatoTheta(unsigned int J, unsigned int nClass, const arma::mat& beta, const arma::mat& Atable);
RcppExport SEXP _ohoegdm_BetatoTheta(SEXP JSEXP, SEXP nClassSEXP, SEXP betaSEXP, SEXP AtableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Atable(AtableSEXP);
    rcpp_result_gen = Rcpp::wrap(BetatoTheta(J, nClass, beta, Atable));
    return rcpp_result_gen;
END_RCPP
}
// computePYaj
arma::mat computePYaj(unsigned int J, unsigned int M, unsigned int nClass, const arma::rowvec& ABETAj, const arma::rowvec& TAUj);
RcppExport SEXP _ohoegdm_computePYaj(SEXP JSEXP, SEXP MSEXP, SEXP nClassSEXP, SEXP ABETAjSEXP, SEXP TAUjSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type ABETAj(ABETAjSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type TAUj(TAUjSEXP);
    rcpp_result_gen = Rcpp::wrap(computePYaj(J, M, nClass, ABETAj, TAUj));
    return rcpp_result_gen;
END_RCPP
}
// computePYa
Rcpp::List computePYa(unsigned int J, unsigned int M, unsigned int nClass, const arma::mat& Atable, const arma::mat& BETA, const arma::mat& TAU);
RcppExport SEXP _ohoegdm_computePYa(SEXP JSEXP, SEXP MSEXP, SEXP nClassSEXP, SEXP AtableSEXP, SEXP BETASEXP, SEXP TAUSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Atable(AtableSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type BETA(BETASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type TAU(TAUSEXP);
    rcpp_result_gen = Rcpp::wrap(computePYa(J, M, nClass, Atable, BETA, TAU));
    return rcpp_result_gen;
END_RCPP
}
// slcm_m2LL
double slcm_m2LL(unsigned int N, unsigned int J, unsigned int M, unsigned int nClass, const arma::mat& Y, const arma::vec& pis, const arma::cube& PY_a);
RcppExport SEXP _ohoegdm_slcm_m2LL(SEXP NSEXP, SEXP JSEXP, SEXP MSEXP, SEXP nClassSEXP, SEXP YSEXP, SEXP pisSEXP, SEXP PY_aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pis(pisSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type PY_a(PY_aSEXP);
    rcpp_result_gen = Rcpp::wrap(slcm_m2LL(N, J, M, nClass, Y, pis, PY_a));
    return rcpp_result_gen;
END_RCPP
}
// Pa1
arma::vec Pa1(unsigned int K, double theta, const arma::vec& lambda0, const arma::vec& lambda1);
RcppExport SEXP _ohoegdm_Pa1(SEXP KSEXP, SEXP thetaSEXP, SEXP lambda0SEXP, SEXP lambda1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lambda0(lambda0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lambda1(lambda1SEXP);
    rcpp_result_gen = Rcpp::wrap(Pa1(K, theta, lambda0, lambda1));
    return rcpp_result_gen;
END_RCPP
}
// slcm_m2LL_HO
double slcm_m2LL_HO(unsigned int N, unsigned int J, unsigned int M, unsigned int nClass, unsigned int K, const arma::mat& Y, const arma::vec& theta, const arma::vec& tau, const arma::vec& lambda, const arma::cube& PY_a, const arma::mat& CLtable);
RcppExport SEXP _ohoegdm_slcm_m2LL_HO(SEXP NSEXP, SEXP JSEXP, SEXP MSEXP, SEXP nClassSEXP, SEXP KSEXP, SEXP YSEXP, SEXP thetaSEXP, SEXP tauSEXP, SEXP lambdaSEXP, SEXP PY_aSEXP, SEXP CLtableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type PY_a(PY_aSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type CLtable(CLtableSEXP);
    rcpp_result_gen = Rcpp::wrap(slcm_m2LL_HO(N, J, M, nClass, K, Y, theta, tau, lambda, PY_a, CLtable));
    return rcpp_result_gen;
END_RCPP
}
// slcm_LLj
double slcm_LLj(unsigned int N, unsigned int M, unsigned int nClass, const arma::vec& Yj, const arma::vec& CLASS, const arma::mat& PY_ajast, const arma::mat& PY_ajtm1);
RcppExport SEXP _ohoegdm_slcm_LLj(SEXP NSEXP, SEXP MSEXP, SEXP nClassSEXP, SEXP YjSEXP, SEXP CLASSSEXP, SEXP PY_ajastSEXP, SEXP PY_ajtm1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Yj(YjSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type CLASS(CLASSSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type PY_ajast(PY_ajastSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type PY_ajtm1(PY_ajtm1SEXP);
    rcpp_result_gen = Rcpp::wrap(slcm_LLj(N, M, nClass, Yj, CLASS, PY_ajast, PY_ajtm1));
    return rcpp_result_gen;
END_RCPP
}
// slcm_LLjm
double slcm_LLjm(unsigned int N, unsigned int M, unsigned int m, unsigned int nClass, const arma::vec& Yj, const arma::vec& CLASS, const arma::mat& PY_ajast, const arma::mat& PY_ajtm1);
RcppExport SEXP _ohoegdm_slcm_LLjm(SEXP NSEXP, SEXP MSEXP, SEXP mSEXP, SEXP nClassSEXP, SEXP YjSEXP, SEXP CLASSSEXP, SEXP PY_ajastSEXP, SEXP PY_ajtm1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type m(mSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Yj(YjSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type CLASS(CLASSSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type PY_ajast(PY_ajastSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type PY_ajtm1(PY_ajtm1SEXP);
    rcpp_result_gen = Rcpp::wrap(slcm_LLjm(N, M, m, nClass, Yj, CLASS, PY_ajast, PY_ajtm1));
    return rcpp_result_gen;
END_RCPP
}
// sampletauast
arma::rowvec sampletauast(unsigned int M, const arma::rowvec& Kaps, double sdMH);
RcppExport SEXP _ohoegdm_sampletauast(SEXP MSEXP, SEXP KapsSEXP, SEXP sdMHSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type Kaps(KapsSEXP);
    Rcpp::traits::input_parameter< double >::type sdMH(sdMHSEXP);
    rcpp_result_gen = Rcpp::wrap(sampletauast(M, Kaps, sdMH));
    return rcpp_result_gen;
END_RCPP
}
// computeLowerBound_Bp
double computeLowerBound_Bp(unsigned int nClass, const arma::mat& LowerAdjTable, const arma::mat& Atable, unsigned int p, const arma::rowvec& betaj, double Bp);
RcppExport SEXP _ohoegdm_computeLowerBound_Bp(SEXP nClassSEXP, SEXP LowerAdjTableSEXP, SEXP AtableSEXP, SEXP pSEXP, SEXP betajSEXP, SEXP BpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type LowerAdjTable(LowerAdjTableSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Atable(AtableSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type betaj(betajSEXP);
    Rcpp::traits::input_parameter< double >::type Bp(BpSEXP);
    rcpp_result_gen = Rcpp::wrap(computeLowerBound_Bp(nClass, LowerAdjTable, Atable, p, betaj, Bp));
    return rcpp_result_gen;
END_RCPP
}
// computeLB
double computeLB(unsigned int nClass, const arma::mat& LBtable, const arma::mat& Atable, unsigned int p, const arma::rowvec& betaj, double Bp, const arma::uvec& Bindices);
RcppExport SEXP _ohoegdm_computeLB(SEXP nClassSEXP, SEXP LBtableSEXP, SEXP AtableSEXP, SEXP pSEXP, SEXP betajSEXP, SEXP BpSEXP, SEXP BindicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type LBtable(LBtableSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Atable(AtableSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type betaj(betajSEXP);
    Rcpp::traits::input_parameter< double >::type Bp(BpSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type Bindices(BindicesSEXP);
    rcpp_result_gen = Rcpp::wrap(computeLB(nClass, LBtable, Atable, p, betaj, Bp, Bindices));
    return rcpp_result_gen;
END_RCPP
}
// identify_check
double identify_check(const arma::mat Q);
RcppExport SEXP _ohoegdm_identify_check(SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(identify_check(Q));
    return rcpp_result_gen;
END_RCPP
}
// lambda_sample
void lambda_sample(unsigned int K, arma::vec& lambda, const arma::vec& m_lam, double v_lam);
RcppExport SEXP _ohoegdm_lambda_sample(SEXP KSEXP, SEXP lambdaSEXP, SEXP m_lamSEXP, SEXP v_lamSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type m_lam(m_lamSEXP);
    Rcpp::traits::input_parameter< double >::type v_lam(v_lamSEXP);
    lambda_sample(K, lambda, m_lam, v_lam);
    return R_NilValue;
END_RCPP
}
// sampletauYast
void sampletauYast(unsigned int N, unsigned int J, unsigned int M, unsigned int nClass, const arma::mat& Y, arma::mat& TAU, arma::mat& Yast, const arma::mat& ABETA, arma::cube& PY_a, const arma::vec& CLASS, arma::vec& MHaccept, double sdMH);
RcppExport SEXP _ohoegdm_sampletauYast(SEXP NSEXP, SEXP JSEXP, SEXP MSEXP, SEXP nClassSEXP, SEXP YSEXP, SEXP TAUSEXP, SEXP YastSEXP, SEXP ABETASEXP, SEXP PY_aSEXP, SEXP CLASSSEXP, SEXP MHacceptSEXP, SEXP sdMHSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type TAU(TAUSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Yast(YastSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ABETA(ABETASEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type PY_a(PY_aSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type CLASS(CLASSSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type MHaccept(MHacceptSEXP);
    Rcpp::traits::input_parameter< double >::type sdMH(sdMHSEXP);
    sampletauYast(N, J, M, nClass, Y, TAU, Yast, ABETA, PY_a, CLASS, MHaccept, sdMH);
    return R_NilValue;
END_RCPP
}
// Q_prime_matrix
arma::mat Q_prime_matrix(unsigned int K, const arma::mat& Atable, const arma::vec& vv);
RcppExport SEXP _ohoegdm_Q_prime_matrix(SEXP KSEXP, SEXP AtableSEXP, SEXP vvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Atable(AtableSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type vv(vvSEXP);
    rcpp_result_gen = Rcpp::wrap(Q_prime_matrix(K, Atable, vv));
    return rcpp_result_gen;
END_RCPP
}
// eta_dina_matrix
arma::mat eta_dina_matrix(const arma::mat& Q);
RcppExport SEXP _ohoegdm_eta_dina_matrix(SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(eta_dina_matrix(Q));
    return rcpp_result_gen;
END_RCPP
}
// pnorm_ln_upper_tail
double pnorm_ln_upper_tail(double& B_p_lowerbound, double& sigma_var_jp);
RcppExport SEXP _ohoegdm_pnorm_ln_upper_tail(SEXP B_p_lowerboundSEXP, SEXP sigma_var_jpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double& >::type B_p_lowerbound(B_p_lowerboundSEXP);
    Rcpp::traits::input_parameter< double& >::type sigma_var_jp(sigma_var_jpSEXP);
    rcpp_result_gen = Rcpp::wrap(pnorm_ln_upper_tail(B_p_lowerbound, sigma_var_jp));
    return rcpp_result_gen;
END_RCPP
}
// q_to_delta
arma::mat q_to_delta(const arma::mat& Q, const arma::mat& Q_prime, unsigned int M);
RcppExport SEXP _ohoegdm_q_to_delta(SEXP QSEXP, SEXP Q_primeSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Q_prime(Q_primeSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(q_to_delta(Q, Q_prime, M));
    return rcpp_result_gen;
END_RCPP
}
// update_slipping_guessing
void update_slipping_guessing(double& slipping, double& guessing, const arma::mat& ab_tilde);
RcppExport SEXP _ohoegdm_update_slipping_guessing(SEXP slippingSEXP, SEXP guessingSEXP, SEXP ab_tildeSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double& >::type slipping(slippingSEXP);
    Rcpp::traits::input_parameter< double& >::type guessing(guessingSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ab_tilde(ab_tildeSEXP);
    update_slipping_guessing(slipping, guessing, ab_tilde);
    return R_NilValue;
END_RCPP
}
// parm_update_nomiss
double parm_update_nomiss(unsigned int N, unsigned int J, unsigned int K, unsigned int nClass, unsigned int M, const arma::mat& Y, arma::mat& Yast, arma::mat& BETA, arma::mat& TAU, arma::vec& CLASS, arma::vec& theta, arma::vec& lambda, arma::vec& tau, arma::mat& Q, arma::mat& DELTA, const arma::mat& Q_prime, const arma::mat& ETA_prime, double& slipping, double& guessing, double omega, const arma::vec& vv, const arma::mat& CLtable, const arma::mat& Atable, const arma::mat& LBtable, const arma::uvec& Bindices, const arma::mat& qtable, unsigned int P, const arma::vec& l1, double m0, const arma::vec& l0, double bq, arma::mat& ABETA, arma::vec& ABETA_sqnorm, arma::cube& PY_a, arma::vec& MHaccept, double sdMH, double& loglike);
RcppExport SEXP _ohoegdm_parm_update_nomiss(SEXP NSEXP, SEXP JSEXP, SEXP KSEXP, SEXP nClassSEXP, SEXP MSEXP, SEXP YSEXP, SEXP YastSEXP, SEXP BETASEXP, SEXP TAUSEXP, SEXP CLASSSEXP, SEXP thetaSEXP, SEXP lambdaSEXP, SEXP tauSEXP, SEXP QSEXP, SEXP DELTASEXP, SEXP Q_primeSEXP, SEXP ETA_primeSEXP, SEXP slippingSEXP, SEXP guessingSEXP, SEXP omegaSEXP, SEXP vvSEXP, SEXP CLtableSEXP, SEXP AtableSEXP, SEXP LBtableSEXP, SEXP BindicesSEXP, SEXP qtableSEXP, SEXP PSEXP, SEXP l1SEXP, SEXP m0SEXP, SEXP l0SEXP, SEXP bqSEXP, SEXP ABETASEXP, SEXP ABETA_sqnormSEXP, SEXP PY_aSEXP, SEXP MHacceptSEXP, SEXP sdMHSEXP, SEXP loglikeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Yast(YastSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type BETA(BETASEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type TAU(TAUSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type CLASS(CLASSSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type DELTA(DELTASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Q_prime(Q_primeSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ETA_prime(ETA_primeSEXP);
    Rcpp::traits::input_parameter< double& >::type slipping(slippingSEXP);
    Rcpp::traits::input_parameter< double& >::type guessing(guessingSEXP);
    Rcpp::traits::input_parameter< double >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type vv(vvSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type CLtable(CLtableSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Atable(AtableSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type LBtable(LBtableSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type Bindices(BindicesSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type qtable(qtableSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type l1(l1SEXP);
    Rcpp::traits::input_parameter< double >::type m0(m0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type l0(l0SEXP);
    Rcpp::traits::input_parameter< double >::type bq(bqSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ABETA(ABETASEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type ABETA_sqnorm(ABETA_sqnormSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type PY_a(PY_aSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type MHaccept(MHacceptSEXP);
    Rcpp::traits::input_parameter< double >::type sdMH(sdMHSEXP);
    Rcpp::traits::input_parameter< double& >::type loglike(loglikeSEXP);
    rcpp_result_gen = Rcpp::wrap(parm_update_nomiss(N, J, K, nClass, M, Y, Yast, BETA, TAU, CLASS, theta, lambda, tau, Q, DELTA, Q_prime, ETA_prime, slipping, guessing, omega, vv, CLtable, Atable, LBtable, Bindices, qtable, P, l1, m0, l0, bq, ABETA, ABETA_sqnorm, PY_a, MHaccept, sdMH, loglike));
    return rcpp_result_gen;
END_RCPP
}
// kappa_initialize
arma::mat kappa_initialize(unsigned int M, unsigned int J);
RcppExport SEXP _ohoegdm_kappa_initialize(SEXP MSEXP, SEXP JSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    rcpp_result_gen = Rcpp::wrap(kappa_initialize(M, J));
    return rcpp_result_gen;
END_RCPP
}
// GenerateAtable
Rcpp::List GenerateAtable(unsigned int nClass, unsigned int K, unsigned int M, unsigned int order);
RcppExport SEXP _ohoegdm_GenerateAtable(SEXP nClassSEXP, SEXP KSEXP, SEXP MSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(GenerateAtable(nClass, K, M, order));
    return rcpp_result_gen;
END_RCPP
}
// QfromD
arma::mat QfromD(unsigned int J, unsigned int K, const arma::mat& DELTA, const arma::mat& DtoQtable);
RcppExport SEXP _ohoegdm_QfromD(SEXP JSEXP, SEXP KSEXP, SEXP DELTASEXP, SEXP DtoQtableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type J(JSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type DELTA(DELTASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type DtoQtable(DtoQtableSEXP);
    rcpp_result_gen = Rcpp::wrap(QfromD(J, K, DELTA, DtoQtable));
    return rcpp_result_gen;
END_RCPP
}
// permuteAtableIndices
arma::vec permuteAtableIndices(unsigned int nClass, unsigned int K, unsigned int M, unsigned int order, const arma::vec& vv, const arma::vec& perm);
RcppExport SEXP _ohoegdm_permuteAtableIndices(SEXP nClassSEXP, SEXP KSEXP, SEXP MSEXP, SEXP orderSEXP, SEXP vvSEXP, SEXP permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type nClass(nClassSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type vv(vvSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type perm(permSEXP);
    rcpp_result_gen = Rcpp::wrap(permuteAtableIndices(nClass, K, M, order, vv, perm));
    return rcpp_result_gen;
END_RCPP
}
// compute_srmr
double compute_srmr(const arma::rowvec& obs_mean, const arma::mat& obs_cov, const arma::rowvec& est_mean, const arma::mat& est_cov);
RcppExport SEXP _ohoegdm_compute_srmr(SEXP obs_meanSEXP, SEXP obs_covSEXP, SEXP est_meanSEXP, SEXP est_covSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type obs_mean(obs_meanSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type obs_cov(obs_covSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type est_mean(est_meanSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type est_cov(est_covSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_srmr(obs_mean, obs_cov, est_mean, est_cov));
    return rcpp_result_gen;
END_RCPP
}
// ohoegdm_cpp
Rcpp::List ohoegdm_cpp(const arma::mat& Y, unsigned int K, unsigned int M, unsigned int order, const arma::vec& l0, const arma::vec& l1, double m0, double bq, double sdMH, unsigned int burnin, unsigned int chain_length);
RcppExport SEXP _ohoegdm_ohoegdm_cpp(SEXP YSEXP, SEXP KSEXP, SEXP MSEXP, SEXP orderSEXP, SEXP l0SEXP, SEXP l1SEXP, SEXP m0SEXP, SEXP bqSEXP, SEXP sdMHSEXP, SEXP burninSEXP, SEXP chain_lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type M(MSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type l0(l0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type l1(l1SEXP);
    Rcpp::traits::input_parameter< double >::type m0(m0SEXP);
    Rcpp::traits::input_parameter< double >::type bq(bqSEXP);
    Rcpp::traits::input_parameter< double >::type sdMH(sdMHSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type chain_length(chain_lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(ohoegdm_cpp(Y, K, M, order, l0, l1, m0, bq, sdMH, burnin, chain_length));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ohoegdm_gen_bijectionvector", (DL_FUNC) &_ohoegdm_gen_bijectionvector, 2},
    {"_ohoegdm_inv_gen_bijectionvector", (DL_FUNC) &_ohoegdm_inv_gen_bijectionvector, 3},
    {"_ohoegdm_CL_gen_invbijection_table", (DL_FUNC) &_ohoegdm_CL_gen_invbijection_table, 3},
    {"_ohoegdm_ETAmat", (DL_FUNC) &_ohoegdm_ETAmat, 4},
    {"_ohoegdm_rTruncNorm_lb", (DL_FUNC) &_ohoegdm_rTruncNorm_lb, 3},
    {"_ohoegdm_rTruncNorm", (DL_FUNC) &_ohoegdm_rTruncNorm, 4},
    {"_ohoegdm_random_Q", (DL_FUNC) &_ohoegdm_random_Q, 2},
    {"_ohoegdm_simSLCM", (DL_FUNC) &_ohoegdm_simSLCM, 8},
    {"_ohoegdm_BetatoTheta", (DL_FUNC) &_ohoegdm_BetatoTheta, 4},
    {"_ohoegdm_computePYaj", (DL_FUNC) &_ohoegdm_computePYaj, 5},
    {"_ohoegdm_computePYa", (DL_FUNC) &_ohoegdm_computePYa, 6},
    {"_ohoegdm_slcm_m2LL", (DL_FUNC) &_ohoegdm_slcm_m2LL, 7},
    {"_ohoegdm_Pa1", (DL_FUNC) &_ohoegdm_Pa1, 4},
    {"_ohoegdm_slcm_m2LL_HO", (DL_FUNC) &_ohoegdm_slcm_m2LL_HO, 11},
    {"_ohoegdm_slcm_LLj", (DL_FUNC) &_ohoegdm_slcm_LLj, 7},
    {"_ohoegdm_slcm_LLjm", (DL_FUNC) &_ohoegdm_slcm_LLjm, 8},
    {"_ohoegdm_sampletauast", (DL_FUNC) &_ohoegdm_sampletauast, 3},
    {"_ohoegdm_computeLowerBound_Bp", (DL_FUNC) &_ohoegdm_computeLowerBound_Bp, 6},
    {"_ohoegdm_computeLB", (DL_FUNC) &_ohoegdm_computeLB, 7},
    {"_ohoegdm_identify_check", (DL_FUNC) &_ohoegdm_identify_check, 1},
    {"_ohoegdm_lambda_sample", (DL_FUNC) &_ohoegdm_lambda_sample, 4},
    {"_ohoegdm_sampletauYast", (DL_FUNC) &_ohoegdm_sampletauYast, 12},
    {"_ohoegdm_Q_prime_matrix", (DL_FUNC) &_ohoegdm_Q_prime_matrix, 3},
    {"_ohoegdm_eta_dina_matrix", (DL_FUNC) &_ohoegdm_eta_dina_matrix, 1},
    {"_ohoegdm_pnorm_ln_upper_tail", (DL_FUNC) &_ohoegdm_pnorm_ln_upper_tail, 2},
    {"_ohoegdm_q_to_delta", (DL_FUNC) &_ohoegdm_q_to_delta, 3},
    {"_ohoegdm_update_slipping_guessing", (DL_FUNC) &_ohoegdm_update_slipping_guessing, 3},
    {"_ohoegdm_parm_update_nomiss", (DL_FUNC) &_ohoegdm_parm_update_nomiss, 37},
    {"_ohoegdm_kappa_initialize", (DL_FUNC) &_ohoegdm_kappa_initialize, 2},
    {"_ohoegdm_GenerateAtable", (DL_FUNC) &_ohoegdm_GenerateAtable, 4},
    {"_ohoegdm_QfromD", (DL_FUNC) &_ohoegdm_QfromD, 4},
    {"_ohoegdm_permuteAtableIndices", (DL_FUNC) &_ohoegdm_permuteAtableIndices, 6},
    {"_ohoegdm_compute_srmr", (DL_FUNC) &_ohoegdm_compute_srmr, 4},
    {"_ohoegdm_ohoegdm_cpp", (DL_FUNC) &_ohoegdm_ohoegdm_cpp, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_ohoegdm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Attributes are assumed binary
// KAPPA is a J by M-2 matrix
// Mixed type items are allowed


//' Generate a vector to map polytomous vector to integers
//'
//' Converts class into a bijection to integers
//'
//' @param K      Number of Attributes
//' @param M      Number of Response Categories
//'
//' @return
//'
//' Return a \eqn{K}-length vector containing the bijection vector.
//'
//' @export
// [[Rcpp::export]]
arma::vec gen_bijectionvector(unsigned int K, unsigned int M)
{
    arma::vec vv(K);
    for (unsigned int k = 0; k < K; k++) {
        vv(k) = pow(M, K - k - 1);
    }
    return vv;
}

// [[Rcpp::export]]
arma::vec inv_gen_bijectionvector(unsigned int K, unsigned int M, double CL)
{
    arma::vec alpha(K);
    for (unsigned int k = 0; k < K; k++) {
        double Mpow = pow(M, K - k - 1);
        double ak = 0.;
        while (((ak + 1.) * Mpow <= CL) & (ak < M)) {
            ak += 1.;
        }
        alpha(k) = ak;
        CL = CL - Mpow * alpha(k);
    }
    return alpha;
}

// function to create a K by nClass table of attribute vectors
// [[Rcpp::export]]
arma::mat CL_gen_invbijection_table(unsigned int K, unsigned int M,
                                    unsigned int nClass)
{
    arma::mat CLtable(K, nClass);
    for (unsigned int cc = 0; cc < nClass; cc++) {
        CLtable.col(cc) = inv_gen_bijectionvector(K, M, cc);
    }
    return CLtable;
}

// Mapping of Q to Delta

// [[Rcpp::export]]
arma::mat ETAmat(unsigned int K, unsigned int J, unsigned int M,

                 const arma::mat &Q)
{

    double nClass = pow(2, K);

    arma::mat ETA(J, nClass);

    for (unsigned int cc = 0; cc < nClass; cc++) {

        arma::vec alpha_c = inv_gen_bijectionvector(K, M, cc);

        for (unsigned int j = 0; j < J; j++) {

            arma::rowvec qj = Q.row(j);
            arma::uvec elementcompare = 1. * (qj.t() <= alpha_c);
            double compare = arma::prod(elementcompare);

            ETA(j, cc) = compare;
        }
    }

    return ETA;
}

//' @title Generate Multinomial Random Variable
//' @description Sample a multinomial random variable for given probabilities.
//' @usage rmultinomial(ps)
//' @param ps A \code{vector} for the probability of each category.
//' @return A \code{vector} from a multinomial with probability ps.
//' @author Steven Andrew Culpepper
//' @noRd
double rmultinomial(const arma::vec &ps, unsigned int M)
{
    double u = R::runif(0, 1);
    arma::vec cps = cumsum(ps);
    arma::vec Ips = arma::zeros<arma::vec>(M);
    Ips.elem(arma::find(cps < u)).fill(1.0);

    return sum(Ips);
}

// [[Rcpp::export]]
double rTruncNorm_lb(double mean, double sd, double b_lb)
{
    double p0 = R::pnorm(b_lb, mean, sd, 1, 0);
    double p1 = 1 - p0;
    double uZ = R::runif(0, 1);
    double Z = R::qnorm(p0 + uZ * p1, mean, sd, 1, 0);
    return (Z);
}

//' @title Generate Dirichlet Random Variable
//' @description Sample a Dirichlet random variable.
//' @usage rDirichlet(deltas)
//' @param deltas A \code{vector} of Dirichlet parameters.
//' @return A \code{vector} from a Dirichlet.
//' @author Steven Andrew Culpepper
//' @noRd
arma::vec rDirichlet(const arma::vec &deltas)
{
    unsigned int C = deltas.n_elem;
    arma::vec Xgamma(C);
    // generating gamma(deltac,1)
    for (unsigned int c = 0; c < C; c++) {
        Xgamma(c) = R::rgamma(deltas(c), 1.0);
    }
    return Xgamma / sum(Xgamma);
}

// [[Rcpp::export]]
double rTruncNorm(double mean, double sd, double w, const arma::vec &ps)
{
    double uZ = R::runif(0, 1);
    double p0 = ps(w);
    double p1 = ps(w + 1);
    double pz = p0 + uZ * (p1 - p0);
    double Z = R::qnorm(pz, mean, sd, 1, 0);
    return Z;
}

// [[Rcpp::export]]
arma::mat random_Q(unsigned int J, unsigned int K)
{
    unsigned int nClass = pow(2, K);
    arma::vec vv = gen_bijectionvector(K, 2);
    arma::vec Q_biject(J);
    Q_biject(arma::span(0, K - 1)) = vv;
    Q_biject(arma::span(K, 2 * K - 1)) = vv;
    arma::vec Jm2K =
        arma::randi<arma::vec>(J - 2 * K, arma::distr_param(1, nClass - 1));
    Q_biject(arma::span(2 * K, J - 1)) = Jm2K;
    Q_biject = arma::shuffle(Q_biject);
    arma::mat Q(J, K);
    for (unsigned int j = 0; j < J; j++) {
        arma::vec qj = inv_gen_bijectionvector(K, 2, Q_biject(j));
        Q.row(j) = qj.t();
    }
    return Q;
}

// M will max Ms when mixed item types are included
//' Simulate Ordinal Item Data from a Sparse Latent Class Model
//' 
//' @param N      Number of Observations
//' @param J      Number of Items
//' @param M      Number of Item Categories (2, 3,  ..., M)
//' @param nClass Number of Latent Classes 
//' @param CLASS  A vector of \eqn{N} observations containing the class ID of the
//'               subject.
//' @param Atable A matrix of dimensions \eqn{M^K \times M^O} containing 
//'               the attribute classes in bijection-form. Note, \eqn{O} refers
//'               to the model's highest interaction order.
//' @param BETA   A matrix of dimensions \eqn{J \times M^K} containing the 
//'               coefficients of the reparameterized \eqn{\beta} matrix.
//' @param KAPPA  A matrix of dimensions \eqn{J \times M} containing the 
//'               category threshold parameters
//' 
//' @return 
//' An ordinal item matrix of dimensions \eqn{N \times J}{N x J} with \eqn{M}
//' response levels.
//' 
//' @seealso [ohoegdm]
//' @export
// [[Rcpp::export]]
arma::mat sim_slcm(unsigned int N, unsigned int J, unsigned int M,
                  unsigned int nClass, const arma::vec &CLASS,
                  const arma::mat &Atable, const arma::mat &BETA,
                  const arma::mat &KAPPA)
{
    arma::cube PY_a = arma::ones<arma::cube>(J, nClass, M + 1);
    PY_a.slice(0) = arma::zeros<arma::mat>(J, nClass);
    // PY_a.slice(M)=arma::ones<arma::mat>(J,nClass);
    for (unsigned int cc = 0; cc < nClass; cc++) {
        arma::rowvec a_alpha = Atable.row(cc);
        for (unsigned int j = 0; j < J; j++) {
            double aBj = arma::accu(a_alpha % BETA.row(j));
            for (unsigned int m = 0; m < M - 1; m++) {
                PY_a(j, cc, m + 1) = R::pnorm(KAPPA(j, m), aBj, 1., 1, 0);
            }
        }
    }
    arma::mat Y(N, J);
    for (unsigned int i = 0; i < N; i++) {
        double class_i = CLASS(i);
        for (unsigned int j = 0; j < J; j++) {
            arma::vec cumulativepsij = PY_a.tube(j, class_i);
            arma::vec psij(M); // not a problem for mixed items; the prob will
                               // be 0 for missing levels
            for (unsigned int m = 0; m < M; m++) {
                psij(m) = cumulativepsij(m + 1) - cumulativepsij(m);
            }
            Y(i, j) = rmultinomial(psij, M);
        }
    }
    return Y;
}

// [[Rcpp::export]]
arma::mat BetatoTheta(unsigned int J, unsigned int nClass,
                      const arma::mat &beta, const arma::mat &Atable)
{
    arma::mat BAp = beta * Atable.t();
    arma::mat theta(J, nClass);
    for (unsigned int j = 0; j < J; j++) {
        for (unsigned int cc = 0; cc < nClass; cc++) {
            theta(j, cc) = R::pnorm(BAp(j, cc), .0, 1., 1, 0);
        }
    }
    return theta;
}

// needs to deal with different Mj
// [[Rcpp::export]]
arma::mat computePYaj(unsigned int J, unsigned int M, unsigned int nClass,
                      const arma::rowvec &ABETAj, const arma::rowvec &KAPPAj)
{

    arma::mat PY_a = arma::ones<arma::mat>(nClass, M + 1);
    PY_a.col(0) = arma::zeros<arma::vec>(nClass);
    // PY_a.col(M)=arma::ones<arma::vec>(nClass);
    // cumulative probs Y given alpha
    for (unsigned int cc = 0; cc < nClass; cc++) {
        double aBj = ABETAj(cc);
        for (unsigned int m = 0; m < M - 1; m++) {
            PY_a(cc, m + 1) = R::pnorm(KAPPAj(m), aBj, 1., 1, 0);
        }
    }
    return PY_a;
}

// need to update this for mixed type items
// probably can delete ABETA and ABETA_sqnorm
// [[Rcpp::export]]
Rcpp::List computePYa(unsigned int J, unsigned int M, unsigned int nClass,
                      const arma::mat &Atable, const arma::mat &BETA,
                      const arma::mat &KAPPA)
{

    arma::cube PY_a = arma::ones<arma::cube>(J, nClass, M + 1);
    PY_a.slice(0) = arma::zeros<arma::mat>(J, nClass);
    // PY_a.slice(M)=arma::ones<arma::mat>(J,nClass);
    arma::mat ABETA(J, nClass);
    arma::vec ABETA_sqnorm = arma::zeros<arma::vec>(nClass);
    // cumulative probs Y given alpha
    // ABETA & ABETA_sqnorm are needed to update alphas
    for (unsigned int cc = 0; cc < nClass; cc++) {
        arma::rowvec a_alpha = Atable.row(cc);
        for (unsigned int j = 0; j < J; j++) {
            double aBj = arma::accu(a_alpha % BETA.row(j));
            ABETA(j, cc) = aBj;
            ABETA_sqnorm(cc) += aBj * aBj;
            for (unsigned int m = 0; m < M - 1; m++) {
                PY_a(j, cc, m + 1) = R::pnorm(KAPPA(j, m), aBj, 1., 1, 0);
            }
        }
    }

    return Rcpp::List::create(Rcpp::Named("ABETA", ABETA),
                              Rcpp::Named("ABETA_sqnorm", ABETA_sqnorm),
                              Rcpp::Named("PY_a", PY_a));
}

// [[Rcpp::export]]
double slcm_m2LL(unsigned int N, unsigned int J, unsigned int M,
                 unsigned int nClass, const arma::mat &Y, const arma::vec &pis,
                 const arma::cube &PY_a)
{

    double m2ll = 0.;
    for (unsigned int i = 0; i < N; i++) {
        arma::rowvec Yi = Y.row(i);
        double Li = 0.;
        for (unsigned int cc = 0; cc < nClass; cc++) {
            double py_a = 1.;
            arma::mat PY_acc = PY_a.subcube(0, cc, 0, J - 1, cc, M);
            for (unsigned int j = 0; j < J; j++) {
                double Yij = Yi(j);
                py_a *= (PY_acc(j, Yij + 1) - PY_acc(j, Yij));
            }
            Li += py_a * pis(cc);
        }
        m2ll += log(Li);
    }
    return -2. * m2ll;
}

// [[Rcpp::export]]
arma::vec Pa1(unsigned int K, double theta, const arma::vec &lambda0,
              const arma::vec &lambda1)
{

    arma::vec ps(K);
    for (unsigned k = 0; k < K; k++) {
        ps(k) = R::pnorm(lambda0(k) + lambda1(k) * theta, 0.0, 1., 1, 0);
    }

    return ps;
}

// [[Rcpp::export]]
double slcm_m2LL_HO(unsigned int N, unsigned int J, unsigned int M,
                    unsigned int nClass, unsigned int K, const arma::mat &Y,
                    const arma::vec &theta, const arma::vec &Tau,
                    const arma::vec &lambda, const arma::cube &PY_a,
                    const arma::mat &CLtable)
{

    double m2ll = 0.;
    for (unsigned int i = 0; i < N; i++) {
        arma::rowvec Yi = Y.row(i);
        double theta_i = theta(i);
        arma::vec Pa1i = Pa1(K, theta_i, Tau, lambda);
        double Li = 0.;
        for (unsigned int cc = 0; cc < nClass; cc++) {
            double py_a = 1.;
            arma::vec ai = CLtable.col(cc);
            arma::mat PY_acc = PY_a.subcube(0, cc, 0, J - 1, cc, M);
            // compute P(alpha'v=c|theta_i)
            double Pacc = 1.;
            for (unsigned int k = 0; k < K; k++) {
                double aik = ai(k);
                double Paik1 = Pa1i(k);
                Pacc *= aik * Paik1 + (1. - aik) * (1. - Paik1);
            }
            for (unsigned int j = 0; j < J; j++) {
                double Yij = Yi(j);
                py_a *= (PY_acc(j, Yij + 1) - PY_acc(j, Yij));
            }
            Li += py_a * Pacc;
        }
        m2ll += log(Li);
    }
    return -2. * m2ll;
}

// [[Rcpp::export]]
double slcm_LLj(unsigned int N, unsigned int M, unsigned int nClass,
                const arma::vec &Yj, const arma::vec &CLASS,
                const arma::mat &PY_ajast, const arma::mat &PY_ajtm1)
{
    double m2ll = 0.;
    for (unsigned int i = 0; i < N; i++) {
        double Yij = Yj(i);
        double class_i = CLASS(i);
        m2ll += log(PY_ajast(class_i, Yij + 1) - PY_ajast(class_i, Yij)) -
                log(PY_ajtm1(class_i, Yij + 1) - PY_ajtm1(class_i, Yij));
    }
    return m2ll;
}

// [[Rcpp::export]]
double slcm_LLjm(unsigned int N, unsigned int M, unsigned int m,
                 unsigned int nClass, const arma::vec &Yj,
                 const arma::vec &CLASS, const arma::mat &PY_ajast,
                 const arma::mat &PY_ajtm1)
{
    double m2ll = 0.;
    for (unsigned int i = 0; i < N; i++) {
        double Yij = Yj(i);
        if ((Yij == m) | (Yij == m + 1)) {
            double class_i = CLASS(i);
            m2ll += log(PY_ajast(class_i, Yij + 1) - PY_ajast(class_i, Yij)) -
                    log(PY_ajtm1(class_i, Yij + 1) - PY_ajtm1(class_i, Yij));
        }
    }
    return m2ll;
}

// [[Rcpp::export]]
arma::rowvec sampleTauast(unsigned int M, const arma::rowvec &Kaps, double sdMH)
{
    arma::vec pk_on(2);
    arma::rowvec KapsMH = Kaps;
    for (unsigned int m = 1; m < M - 1; m++) {
        double uMH = R::runif(0.0, 1.0);
        if (m == 1) {
            // pk_on(0) = 0.;
            pk_on(0) = R::pnorm(KapsMH(m - 1), Kaps(m), sdMH, 1, 0);
            pk_on(1) = R::pnorm(Kaps(m + 1), Kaps(m), sdMH, 1, 0);
        } else if (m == M - 2) {
            pk_on(0) = R::pnorm(KapsMH(m - 1), Kaps(m), sdMH, 1, 0);
            pk_on(1) = 1.;
        } else {
            pk_on(0) = R::pnorm(KapsMH(m - 1), Kaps(m), sdMH, 1, 0);
            pk_on(1) = R::pnorm(Kaps(m + 1), Kaps(m), sdMH, 1, 0);
        }
        KapsMH(m) = R::qnorm(pk_on(0) + uMH * (pk_on(1) - pk_on(0)), Kaps(m),
                             sdMH, 1, 0);
    }
    return KapsMH;
}

// this functions needs to be checked for generalizability to order<nClass
// [[Rcpp::export]]
double computeLowerBound_Bp(unsigned int nClass, const arma::mat &LowerAdjTable,
                            const arma::mat &Atable, unsigned int p,
                            const arma::rowvec &betaj, double Bp)
{
    // compute A'Beta_j
    arma::vec Abetaj = Atable * betaj.t();
    // Extract column p of Atable design matrix
    arma::vec Ap = Atable.col(p);
    // Find the elements that include Bp in the computation of Abetaj
    // by finding the elements in Ap that are equal to one
    arma::uvec elements_equal_one = find(Ap == 1);
    // Count the number of elements equal to one
    unsigned int numberofones = elements_equal_one.n_elem;
    // Create an object to store the lower bounds for each case
    arma::vec lbparts(numberofones);
    // loop over the cases with Bp in the computation of Abetaj to find each max
    // lowerbound
    for (unsigned int r = 0; r < numberofones; r++) {
        // translate the r index to the original index
        double rr = elements_equal_one(r);
        // subtract Bp from Abetaj for case r, arr(p)'betaj(p) = arr'betaj-Bjp
        double Abetajr_minus_beta_jp = Abetaj(rr) - Bp;
        // Flag elements of Abetaj that are adjacent to case r, less than it,
        // and have Ap = 0 these rows are a_adj with latent mean a_adj'betaj
        arma::vec adjacentflag = LowerAdjTable.col(rr) % (1 - Ap);
        // Find the elements
        arma::uvec adjacentindices = find(adjacentflag == 1);
        // compute the lowerbound as max of Abetaj for adjacent cases less
        // Abetajr_minus_beta_jp
        // Bjp > max a_adj'betaj - arr(p)'betaj(p)
        lbparts(r) = arma::max(Abetaj(adjacentindices)) - Abetajr_minus_beta_jp;
    }

    return arma::max(lbparts);
}

// [[Rcpp::export]]
double computeLB(unsigned int nClass, const arma::mat &LBtable,
                 const arma::mat &Atable, unsigned int p,
                 const arma::rowvec &betaj, double Bp,
                 const arma::uvec &Bindices)
{
    double Bplb = -100.;
    if (p < nClass - 1) {
        arma::vec lbmax(2);
        arma::vec LBp = LBtable.col(p);
        arma::uvec lbinds = find(LBp == 1); // find lower classes
        arma::rowvec ap =
            Atable.row(Bindices(p)); // find ap design vector for predictor p
        arma::rowvec betajpeq0 = betaj;
        betajpeq0(p) = 0.; // create a vector of betaj with element p = 0
        double gamp =
            arma::dot(ap, betajpeq0); // compute gammap without predictor p
        arma::vec gams = (Atable.rows(lbinds)) *
                         betaj.t(); // compute gammas for all lower classes
        lbmax(0) = arma::max(gams) - gamp; // compute diff bw max of lower
                                           // classes and gammap without betap
        arma::uvec lbinds2 = find(LBp == 2); // find greater or equal classes
        arma::vec gamdiff = gamp - (Atable.rows(lbinds2)) * betajpeq0.t();
        lbmax(1) = arma::max(gamdiff);
        Bplb = arma::max(lbmax);
    }
    if (p == nClass - 1) { // need this because there is no 2 in LBtable for
                           // last coefficient
        arma::vec LBp = LBtable.col(p);
        arma::uvec lbinds = find(LBp == 1); // find lower classes
        arma::rowvec ap = Atable.row(Bindices(p));
        arma::rowvec betajpeq0 = betaj;
        betajpeq0(p) = 0.;
        double gamp = arma::accu(ap % betajpeq0);
        arma::vec gams = (Atable.rows(lbinds)) * betaj.t();
        Bplb = arma::max(gams) - gamp;
    }
    return Bplb;
}

// [[Rcpp::export]]
double identify_check(const arma::mat Q)
{
    unsigned int K = Q.n_cols;
    unsigned int J = Q.n_rows;

    arma::mat ones_zero_on_diag = -1 * arma::ones<arma::mat>(K, K);
    arma::vec zeros_K = arma::zeros<arma::vec>(K);
    ones_zero_on_diag.diag() = zeros_K;
    arma::vec c_sum = (arma::sum(Q, 0)).t();
    arma::vec r_sum = arma::sum(Q, 1);
    arma::mat I_check = Q * ones_zero_on_diag;
    arma::mat I_count = arma::zeros<arma::mat>(J, K);
    I_count.elem(arma::find(I_check > -1)).fill(1.0);
    arma::vec n_ek = (arma::sum(I_count, 0)).t();

    double min_c = (arma::min(c_sum) > 2);
    double min_r = (arma::min(r_sum) > 0);
    double min_ek = (arma::min(n_ek) > 1);

    return (min_c + min_r + min_ek > 2);
}

// note this freely samples lambda, so attribute levels may flip, but the
// estimation of B may be better
// [[Rcpp::export]]
void lambda_sample(unsigned int K, arma::vec &lambda, const arma::vec &m_lam,
                   double v_lam)
{
    double sd_lam = sqrt(v_lam);
    for (unsigned int k = 0; k < K; k++) {
        lambda(k) = R::rnorm(m_lam(k), sd_lam);
    }
}

// // [[Rcpp::export]]
// void  lambda_sample(unsigned int K,arma::vec& lambda,const arma::vec&
// m_lam,double v_lam){
//   double sd_lam =sqrt(v_lam);
//   for(unsigned int k=0;k<K;k++){
//     double uk = R::runif(0,1);
//     double p0 = R::pnorm(0.,m_lam(k),sd_lam,1,0);
//     double p1 = 1.-p0;
//     lambda(k) = R::qnorm(p0+uk*p1,m_lam(k),sd_lam,1,0);
//   }
// }

// [[Rcpp::export]]
void sampleTauYast(unsigned int N, unsigned int J, unsigned int M,
                   unsigned int nClass, const arma::mat &Y, arma::mat &KAPPA,
                   arma::mat &Yast, const arma::mat &ABETA, arma::cube &PY_a,
                   const arma::vec &CLASS, arma::vec &MHaccept, double sdMH)
{
    arma::vec pk_on(4);
    arma::vec pky_on(4);
    for (unsigned int j = 0; j < J; j++) {
        arma::vec Yj = Y.col(j);
        arma::rowvec ABETAj = ABETA.row(j);
        arma::rowvec Kaps = KAPPA.row(j);
        // Sample MH Threshold candidates
        arma::rowvec KapsMH = sampleTauast(M, Kaps, sdMH);
        // Step 1a: compute m part related to thresholds
        double lnmpart = .0;
        for (unsigned int m = 1; m < M - 1; m++) {
            pk_on(0) = R::pnorm(KapsMH(m - 1), Kaps(m), sdMH, 1, 0);
            pk_on(2) = R::pnorm(Kaps(m - 1), KapsMH(m), sdMH, 1, 0);
            if (m == M - 2) {
                pk_on(1) = 1.;
                pk_on(3) = 1.;
            } else {
                pk_on(1) = R::pnorm(Kaps(m + 1), Kaps(m), sdMH, 1, 0);
                pk_on(3) = R::pnorm(KapsMH(m + 1), KapsMH(m), sdMH, 1, 0);
            }
            lnmpart += log(pk_on(1) - pk_on(0)) - log(pk_on(3) - pk_on(2));
        }
        // Step 1b: compute i part related to individuals
        arma::mat PYajMH = computePYaj(J, M, nClass, ABETAj, KapsMH);
        arma::mat PYaj = PY_a.subcube(j, 0, 0, j, nClass - 1, M);
        double lnipart = .0;
        for (unsigned int i = 0; i < N; i++) {
            double class_i = CLASS(i);
            double Yij = Yj(i);
            pky_on(0) = PYajMH(class_i, Yij);
            pky_on(1) = PYajMH(class_i, Yij + 1.);
            pky_on(2) = PYaj(class_i, Yij);
            pky_on(3) = PYaj(class_i, Yij + 1.);
            lnipart += log(pky_on(1) - pky_on(0)) - log(pky_on(3) - pky_on(2));
        }
        // Step 1c: compute part related to prior, prior var = 9
        // double
        // lnpriorpart=-.5*(arma::dot(KapsMH,KapsMH)-arma::dot(Kaps,Kaps))/9.;
        // Step 2: Compute R and min prob for rejection & Update Kappas and Z
        // NOte that R_MH = mpart * ipart *priorpart
        double lnuR = log(R::runif(0.0, 1.0));
        if (lnuR < lnmpart + lnipart) { //+ lnpriorpart
            MHaccept(j) = 1.;
            Kaps = KapsMH;
            KAPPA.row(j) = Kaps;
            PYaj = PYajMH;
            // arma::mat PYaj = PY_a.subcube(j,0,0,j,nClass-1,M);
            PY_a.subcube(j, 0, 0, j, nClass - 1, M) = PYaj;
            for (unsigned int i = 0; i < N; i++) {
                double class_i = CLASS(i);
                double Yij = Yj(i);
                double aiBj = ABETAj(class_i);
                arma::rowvec pYij = PYaj.row(class_i);
                double Zijk = rTruncNorm(aiBj, 1., Yij, pYij.t());
                Yast(i, j) = Zijk;
            }
        } else {
            MHaccept(j) = .0;
        }
    }
    // return MHaccept;
}




// [[Rcpp::export]]
arma::mat Q_prime_matrix(unsigned int K, const arma::mat &Atable,
                         const arma::vec &vv)
{
    
    // The number of rows for Q prime must always be p
    // P here is found under the number of columns not the rows of the
    // A table
    // Nuance: Model is not saturated, we may need to update the code elsewhere
    arma::mat Q_prime(Atable.n_cols, K);
    
    // Setup a basis vector
    arma::vec e_k = arma::zeros<arma::vec>(K);
    
    for (unsigned int k = 0; k < K; ++k) {
        
        // Fill the k-th position with 1
        e_k(k) = 1;
        
        // Perform a dot product
        unsigned int col_ind = arma::dot(e_k, vv);
        
        // Retrieve from the Atable the appropriate column based off of the
        // matrix bijection vector vv.
        Q_prime.col(k) = Atable.col(col_ind);
        
        // Reset the k-th position for next iteration
        e_k(k) = 0;
    }
    
    return Q_prime;
}

// [[Rcpp::export]]
arma::mat eta_dina_matrix(const arma::mat &Q)
{
    
    // Set up attributes for populating eta.
    unsigned int K = Q.n_cols, J = Q.n_rows;
    
    // Compute the total number of attribute classes C = 2^K
    double nClass = pow(2, K);
    
    // Setup storage for the eta
    arma::mat eta_jc(J, nClass);
    
    for (unsigned int cc = 0; cc < nClass; ++cc) {
        arma::vec alpha_c = inv_gen_bijectionvector(K, 2, cc);
        
        for (unsigned int j = 0; j < J; ++j) {
            arma::rowvec qj = Q.row(j);
            // Switch to as_scalar
            double compare = arma::as_scalar(qj * alpha_c - qj * qj.t());
            eta_jc(j, cc) = (compare >= 0);
        }
    }
    
    return eta_jc;
}

// [[Rcpp::export]]
double pnorm_ln_upper_tail(double &B_p_lowerbound, double &sigma_var_jp) {
    return R::pnorm(B_p_lowerbound / sqrt(sigma_var_jp), 0.0, 1.0, 0, 1);
}

// [[Rcpp::export]]
arma::mat q_to_delta(const arma::mat& Q,
                     const arma::mat& Q_prime,
                     unsigned int M) {
    
    unsigned int K = Q.n_cols;
    
    // Table of deltas from Iteration 1 back.
    arma::mat ETA_prime = eta_dina_matrix(Q_prime);
    
    // Create a bijection vector
    arma::vec vv = gen_bijectionvector(K, M);
    
    // J x 2^K (or lower order with J x 2^P)
    arma::mat DELTA(Q.n_rows, Q_prime.n_rows);
    
    for(unsigned int j = 0; j < Q.n_rows; ++j) {
        // Translate to an integer, pick out column of ETA prime and store it
        // as a row in DELTA
        DELTA.row(j) = ETA_prime.col( dot(Q.row(j), vv) ).t();
    }
    
    return DELTA;
}

// [[Rcpp::export]]
void update_slipping_guessing(double &slipping, double &guessing,
                              const arma::mat &ab_tilde)
{
    // Global Slipping and Guessing update across parameters
    
    // Initialize on a uniform
    double slipping_unif = R::runif(0.0, 1.0);
    double guessing_unif = R::runif(0.0, 1.0);
    
    // Retrieve old slipping value for item j
    double slipping_old = slipping;
    
    // Draw guessing conditioned upon slipping - 1
    double ab_g1 = ab_tilde(0, 1); // eta 0, delta 1
    double ab_g0 = ab_tilde(0, 0); // eta 0, delta 0
    
    double pg = R::pbeta(1.0 - slipping_old, ab_g1 + 1., ab_g0 + 1., 1, 0);
    
    double guessing_new =
        R::qbeta(guessing_unif * pg, ab_g1 + 1., ab_g0 + 1., 1, 0);
    
    // Draw slipping conditioned upon guessing
    double ab_s1 = ab_tilde(1, 1); // eta_prime 1, delta 1
    double ab_s0 = ab_tilde(1, 0); // eta_prime 1, delta 0
    
    double ps = R::pbeta(1.0 - guessing_new, ab_s0 + 1., ab_s1 + 1., 1, 0);
    
    double slipping_new =
        R::qbeta(slipping_unif * ps, ab_s0 + 1., ab_s1 + 1., 1, 0);
    
    // Update slipping and guessing for all items
    slipping = slipping_new;
    guessing = guessing_new;
}


// [[Rcpp::export]]
double parm_update_nomiss(
    unsigned int N, unsigned int J, unsigned int K, unsigned int nClass,
    unsigned int M, const arma::mat &Y, arma::mat &Yast, arma::mat &BETA,
    arma::mat &KAPPA, arma::vec &CLASS, arma::vec &theta, arma::vec &lambda,
    arma::vec &Tau, arma::mat &Q, 
    arma::mat &DELTA,
    const arma::mat &Q_prime,
    const arma::mat &ETA_prime,
    double &slipping, double &guessing,
    double omega, const arma::vec &vv,
    const arma::mat &CLtable, const arma::mat &Atable, const arma::mat &LBtable,
    const arma::uvec &Bindices, const arma::mat &qtable, unsigned int P,
    const arma::vec &l1, double m0, const arma::vec &l0, 
    double bq, arma::mat &ABETA, arma::vec &ABETA_sqnorm, arma::cube &PY_a,
    arma::vec &MHaccept, double sdMH, double &loglike)
{

    // Rows are Eta and Columns Delta
    arma::mat ab_tilde = arma::zeros<arma::mat>(2, 2);
    
    // update KAPPA and Yast
    sampleTauYast(N, J, M, nClass, Y, KAPPA, Yast, ABETA, PY_a, CLASS, MHaccept,
                  sdMH);

    // update classes + store info for Betas and pi
    arma::mat ApA = arma::zeros<arma::mat>(P, P);
    arma::mat ApZ = arma::zeros<arma::mat>(P, J);
    arma::vec Aastmthetalambda = arma::zeros<arma::vec>(K);
    arma::vec Aasttheta = arma::zeros<arma::vec>(K);

    double v_theta = 1. / (arma::dot(lambda, lambda) + 1.);
    double sd_theta = sqrt(v_theta);
    loglike = 0.;

    for (unsigned int i = 0; i < N; i++) {
        arma::rowvec Yi = Y.row(i);
        // double class_i=CLASS(i);
        double theta_i = theta(i);
        arma::vec numerator(nClass);
        double denominator = 0.;
        arma::vec Pa1i = Pa1(K, theta_i, Tau, lambda);
        for (unsigned int cc = 0; cc < nClass; cc++) {
            double picc = 1.;
            arma::vec ai = CLtable.col(cc);
            double Pacc = 1.;
            for (unsigned int k = 0; k < K; k++) {
                double aik = ai(k);
                double Paik1 = Pa1i(k);
                Pacc *= aik * Paik1 + (1. - aik) * (1. - Paik1);
            }
            for (unsigned int j = 0; j < J; j++) {
                double Yij = Yi(j);
                picc *= (PY_a(j, cc, Yij + 1.) - PY_a(j, cc, Yij));
            }
            numerator(cc) = picc * Pacc;
            denominator += picc * Pacc;
        }
        arma::vec pai = numerator / denominator;
        double class_i = rmultinomial(pai, nClass);
        CLASS(i) = class_i;
        arma::rowvec a_alpha = Atable.row(class_i);
        loglike += log(denominator);
        // arma::rowvec Yasti(J);
        //
        // //sample Y augmented data
        // for(unsigned int j=0;j<J;j++){
        //   double Yij = Yi(j);
        //   double aiBj=ABETA(j,class_i);
        //   arma::vec pYij= PY_a.tube(j,class_i);
        //   double Zijk= rTruncNorm(aiBj,1.,Yij,pYij);
        //   Yasti(j) = Zijk;
        // }
        // //Yast.row(i)=Yasti;
        arma::rowvec Yasti = Yast.row(i);
        ApA += a_alpha.t() * a_alpha;
        ApZ += a_alpha.t() * Yasti;

        // sample alpha augmented data
        arma::vec ai = CLtable.col(class_i);
        arma::vec mean_i(K);
        arma::vec aiast(K);
        for (unsigned int k = 0; k < K; k++) {
            double aik = ai(k);
            double paik1 = Pa1i(k); // R::pnorm(0.0,mean_i(k),1.,1,0);
            arma::vec paik(3);
            paik(0) = 0.;
            paik(1) = 1. - paik1;
            paik(2) = 1.;
            mean_i(k) = Tau(k) + lambda(k) * theta_i;
            double aikast = rTruncNorm(mean_i(k), 1., aik, paik);
            aiast(k) = aikast;
        }
        // update theta
        arma::vec Wtilde_lam = aiast - Tau;
        double m_theta = v_theta * arma::dot(Wtilde_lam, lambda);
        theta_i = R::rnorm(m_theta, sd_theta);

        Aastmthetalambda += aiast - lambda * theta_i;
        Aasttheta += aiast * theta_i;
        theta(i) = theta_i;
    }
    // arma::vec oneN=arma::ones<arma::vec>(N);
    double s2Tau_inv = 3.;
    // update Tau + s2Tau_inv
    double v_Tau = 1. / (s2Tau_inv + double(N));
    // arma::rowvec m_Taus = v_Tau*arma::sum( W-theta*lambda.t(),0);
    // Tau = arma::randn<arma::vec>(K)*sqrt(v_Tau)+m_Taus.t();
    arma::vec m_Taus = v_Tau * Aastmthetalambda;
    Tau = arma::randn<arma::vec>(K) * sqrt(v_Tau) + m_Taus;
    // arma::rowvec m_Taus = v_Tau*arma::sum( Aast-theta*lambda.t(),0);
    // Tau = arma::randn<arma::vec>(K)*sqrt(v_Tau)+m_Taus.t();

    // double s2Tau_inv_new = 1.;
    // update lambda + s2lam_inv
    double s2lam_inv = 2;
    double v_lam = 1. / (arma::dot(theta, theta) + s2lam_inv); // /double(K);
    // arma::mat Wtilde_lam = W-oneN*Tau.t();
    // arma::vec m_lams=v_lam*(Wtilde_lam.t()*theta);
    arma::vec AastmTautheta = Aasttheta - arma::accu(theta) * Tau;
    arma::vec m_lams = v_lam * AastmTautheta;
    // arma::vec lambda_k(1);
    // arma::vec oneK = arma::ones<arma::vec>(K);
    // lambda_sample(1,lambda_k,oneK.t()*m_lams,v_lam);
    // lambda = oneK*lambda_k;
    lambda_sample(K, lambda, m_lams, v_lam);

    // update Q, BETA, this version uses MH for q. v2_1 uses Gibbs
    for (unsigned int j = 0; j < J; j++) {
        arma::vec ApZj = ApZ.col(j);
        arma::rowvec betaj = BETA.row(j);
        arma::rowvec q_j = Q.row(j);
        
        // ETA 2^K to P matrix
        // Retrieve: 1 x P
        arma::rowvec delta_j = DELTA.row(j);
        
        for (unsigned int k = 0; k < K; k++) {
            
            // qjk = Q(j,k);
            // checking whether 1 is possible
            
            arma::mat Q1 = Q;
            Q1(j, k) = 1;
            
            double flag1 = identify_check(Q1);
            if (flag1 == 1) {
                
                // checking whether 0 is possible
                arma::mat Q0 = Q;
                Q0(j, k) = 0;
                double flag0 = identify_check(Q0);
                
                // update based upon posterior for qjk
                if (flag0 == 1) {
                    arma::uvec pfork = find(qtable.row(k) == 1);
                    
                    // find qast,qtm1
                    double qtm1 = q_j(k);
                    double qast = 1. - qtm1;
                    
                    // Determine the current Q matrix eta column
                    double col_t_eta_mat = arma::dot(q_j, vv);
                    
                    // Determine the proposed Q matrix eta column
                    q_j(k) = 1 - q_j(k);
                    double col_ast_eta_mat = arma::dot(q_j, vv);
                    q_j(k) = 1 - q_j(k);
                    
                    // Perform counts with the delta matrix ----
                    
                    // Relate the Delta matrix to Q prime
                    arma::vec related_q_prime_delta =
                        delta_j.t() % Q_prime.col(k);
                    
                    arma::vec related_q_prime_1_minus_delta =
                        (1 - delta_j.t()) % Q_prime.col(k);
                    
                    arma::vec ETA_prime_previous = ETA_prime.col(col_t_eta_mat);
                    arma::vec ETA_prime_proposed = ETA_prime.col(col_ast_eta_mat);
                    
                    // Count slipping components ----
                    
                    // Previous Q matrix entry in chain 1 - slipping count
                    double count_1_minus_slipping_previous = arma::dot(
                        related_q_prime_delta, ETA_prime_previous);
                    
                    // Proposed Q matrix entry in chain 1 - slipping count
                    double count_1_minus_slipping_proposed = arma::dot(
                        related_q_prime_delta, ETA_prime_proposed);
                    
                    // Previous Q matrix entry in chain slipping count
                    double count_slipping_previous =
                        arma::dot(related_q_prime_1_minus_delta,
                                  ETA_prime_previous);
                    
                    // Proposed Q matrix entry in chain slipping count
                    double count_slipping_proposed =
                        arma::dot(related_q_prime_1_minus_delta,
                                  ETA_prime_proposed);
                    
                    // Count guessing components ----
                    
                    // Previous Q matrix entry in chain guessing count
                    double count_guessing_previous =
                        arma::dot(related_q_prime_delta,
                                  1.0 - ETA_prime_previous);
                    
                    // Proposed Q matrix entry in chain guessing count
                    double count_guessing_proposed =
                        arma::dot(related_q_prime_delta,
                                  1.0 - ETA_prime_proposed);
                    
                    // Previous Q matrix entry in chain guessing count
                    double count_1_minus_guessing_previous =
                        arma::dot(related_q_prime_1_minus_delta,
                                  1.0 - ETA_prime_previous);
                    
                    // Proposed Q matrix entry in chain guessing count
                    double count_1_minus_guessing_proposed =
                        arma::dot(related_q_prime_1_minus_delta,
                                  1.0 - ETA_prime_proposed);
                    
                    // Perform the Metropolis Hastings Update
                    // Sample q_jk from conditional Bernoulli
                    double u = R::runif(0, 1);
                    
                    // Natural log variant of the acceptance ratio
                    double ln_acceptance_ratio =
                        (count_1_minus_slipping_proposed -
                        count_1_minus_slipping_previous) *
                        log(1.0 - slipping) +
                        (count_guessing_proposed - count_guessing_previous) *
                        log(guessing) +
                        (count_slipping_proposed - count_slipping_previous) *
                        log(slipping) +
                        (count_1_minus_guessing_proposed -
                        count_1_minus_guessing_previous) *
                        log(1.0 - guessing) +
                        (qast - qtm1) * (log(omega) - log(1.0 - omega));
                    
                    double flag =
                        1. * (ln_acceptance_ratio > log(u));       // MH update
                    double qjk = qast * flag + (1. - flag) * qtm1; // MH update
                    
                    Q(j, k) = qjk;
                    q_j(k) = qjk;
                    
                } // end if metropolis hasting update
                
            } // end if check for proposed Q matrix identifiability.
            
        } // end for loop for Q matrix K update
        
        // sumQ+=arma::accu(qj);
        // Q.row(j)=qj;
        
        // Retrieve the ETA prime row based on the updated q vector
        arma::vec ETA_prime_q_j = ETA_prime.col(arma::dot(q_j, vv));
        
        // P is just 0 to 2^(K-1) [saturated model]
        arma::uvec idx = arma::linspace<arma::uvec>(0, P - 1, P);
        
        
        for (unsigned int p = 0; p < P; p++) {
            
            // Update the Delta ----
            
            // Initialize the delta_jp to 1 because the intercept is always
            // held at 1 since we want it in the model.
            double delta_jp = 1;
            
            // Setup variables to compute the lower bound
            // This is only computed when we are not on the intercept update.
            double B_p = 0;
            double B_p_lowerbound = 0;
            
            // Retrieve the p-th value from ETA_prime_j
            double ETA_prime_q_jp = ETA_prime_q_j(p);
            
            // Silently the intercept is going to be 1 automatically, so
            // no update will be performed
            if (p > 0) {
                
                // Obtain delta_jp values
                double delta_jp_previous = delta_j(p);
                double delta_jp_proposed = 1.0 - delta_jp_previous;
                
                //  Obtain the variance
                double sigma_var_jp_previous =
                    delta_jp_previous * 1.0 / l1(p) +
                    (1. - delta_jp_previous) * 1.0 / l0(p);
                
                double sigma_var_jp_proposed =
                    delta_jp_proposed * 1.0 / l1(p) +
                    (1.0 - delta_jp_proposed) * 1.0 / l0(p);
                
                // Compute the lower bound
                B_p = betaj(p);
                B_p_lowerbound =
                    computeLB(nClass, LBtable, Atable, p, betaj, B_p, Bindices);
                
                // Obtain the standard normal probability in a left tail
                // distribution on the log normal scale.
                // Equivalent to ln(1 - psi())
                double ln_normalizing_constant_previous =
                    pnorm_ln_upper_tail(B_p_lowerbound, sigma_var_jp_previous);
                
                double ln_normalizing_constant_proposed =
                    pnorm_ln_upper_tail(B_p_lowerbound, sigma_var_jp_proposed);
                
                // Compute the natural log ratios for both components in delta
                // MH update
                double ln_ratio_betas = 0.5 * (log(sigma_var_jp_previous) -
                                               log(sigma_var_jp_proposed)) +
                                               (ln_normalizing_constant_previous -
                                               ln_normalizing_constant_proposed) -
                                               0.5 *
                                               (1.0 / sigma_var_jp_proposed -
                                               1.0 / sigma_var_jp_previous) *
                                               B_p * B_p;
                
                double ln_ratio_deltas =
                    (delta_jp_proposed - delta_jp_previous) *
                    (ETA_prime_q_jp * log((1.0 - slipping) / slipping) +
                    (1 - ETA_prime_q_jp) * log(guessing / (1.0 - guessing)));
                
                double ln_acceptance_delta_update =
                    ln_ratio_betas + ln_ratio_deltas;
                
                // Perform the Metropolis Hastings Update
                // Sampling Delta_jp from a truncated normal using MH
                double u = R::runif(0, 1);
                
                double flag_delta_update =
                    1. * (ln_acceptance_delta_update > log(u));
                
                // Based on MH update change the delta or keep it as-is
                delta_jp = delta_jp_proposed * flag_delta_update +
                    (1. - flag_delta_update) * delta_jp_previous;
                delta_j(p) = delta_jp;
                DELTA(j, p) = delta_jp;
                
                // Count only values not in the intercept based on:
                // ETA, DELTA
                ab_tilde(ETA_prime_q_jp, delta_jp) += 1;
                
            } //
            
            
            double djp = delta_jp; // deltaj(p);
            double vforb = djp * l1(p) + (1. - djp) * l0(p);
            double ApApp = ApA(p, p);
            double vb = 1. / (ApApp + vforb);
            double ApAb = arma::accu(ApA.row(p) % betaj) - ApApp * betaj(p);
            double mub = vb * (ApZj(p) - ApAb);
            double sdb = sqrt(vb);
            if (p == 0) {
                mub += vb * vforb * m0;
                betaj(p) = R::rnorm(mub, sqrt(vb));
            } else {
                double Bp = betaj(p);
                double Bplb =
                    computeLB(nClass, LBtable, Atable, p, betaj, Bp, Bindices);
                // betaj(p) = rTruncNorm_b(mub,sqrt(vb),1.);
                if ((Bplb - mub) / sdb > 4.264891) {
                    betaj(p) = Bplb;
                } else {
                    betaj(p) = rTruncNorm_lb(mub, sdb, Bplb);
                }
            }
        }
        
        arma::rowvec ABETAj = (Atable * betaj.t()).t();
        arma::rowvec KAPPAj = KAPPA.row(j);
        arma::mat PYajtm1 = computePYaj(J, M, nClass, ABETAj, KAPPAj);
        PY_a.subcube(j, 0, 0, j, nClass - 1, M) = PYajtm1;
        ABETA.row(j) = ABETAj;
        BETA.row(j) = betaj;
    }

    // update ABETA_sqnorm
    // arma::vec ABETAsqtmp=arma::zeros<arma::vec>(nClass);
    // for(unsigned int cc=0;cc<nClass;cc++){
    //   arma::rowvec a_alpha=Atable.row(cc);
    //   for(unsigned int j=0;j<J;j++){
    //     double aBj = arma::accu(a_alpha%BETA.row(j));
    //     ABETAsqtmp(cc)+=aBj*aBj;
    //   }
    // }
    // ABETA_sqnorm=ABETAsqtmp;
    
    // Perform the slipping and guessing update
    update_slipping_guessing(slipping, guessing, ab_tilde);
    
    double sumQ = arma::accu(Q);
    double omeganew = R::rbeta(sumQ + 1.0, double(J * K) - sumQ + bq);
    return omeganew;
}

// [[Rcpp::export]]
arma::mat kappa_initialize(unsigned int M, unsigned int J)
{

    arma::mat KAP0(J, M - 1);
    (KAP0.col(0)).fill(.0);
    if (M > 2) {
        for (unsigned int j = 0; j < J; j++) {
            for (unsigned int m = 1; m < M - 1; m++) {
                KAP0(j, m) = KAP0(j, m - 1) + R::runif(.8, 1.2);
            }
        }
    }
    return KAP0;
}

// Note the M in this function is Malpha, NOT the number of response categories

//' Generate tables that store different design elements
//' 
//' Each table provides a "cache" of pre-computed values.
//' 
//' @param nClass   Number of Attribute Classes
//' @param M        Number of Responses
//' @param K        Number of Attributes
//' @param order    Highest interaction order to consider. 
//'                 Default model-specified `k`.
//'
//' @return
//' Return a `list` containing the table caches for different parameters
//'
//' @details
//' This is **an internal function** briefly used to simulate data and, thus, has
//' been exported into _R_ as well as documented. **Output from this function can
//' change in future versions.**
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List GenerateAtable(unsigned int nClass, unsigned int K, unsigned int M,
                          unsigned int order)
{
    arma::mat FullAtable(nClass, nClass);
    arma::mat FullLBtable(nClass, nClass);
    arma::mat FullDtoQtable(K, nClass);
    arma::mat Fulladjtable(nClass, nClass);
    arma::vec model(nClass);
    for (unsigned int cr = 0; cr < nClass; cr++) {
        arma::vec alpha_r = inv_gen_bijectionvector(K, M, cr);
        double nof0s = 0.;
        for (unsigned int k = 0; k < K; k++) {
            nof0s += 1. * (alpha_r(k) == 0);
            FullDtoQtable(k, cr) = 1. * (alpha_r(k) > 0);
        }
        model(cr) = 1. * (nof0s > double(K - order) - 1);
        for (unsigned int cc = 0; cc < nClass; cc++) {
            arma::vec alpha_c = inv_gen_bijectionvector(K, M, cc);
            double mindiff = arma::min(alpha_r - alpha_c);
            FullAtable(cr, cc) = 1. * (mindiff > -1);
            double maxdiff = arma::accu(abs(alpha_c - alpha_r));
            FullLBtable(cr, cc) = 1. * (maxdiff == 1) * (mindiff < 0) +
                                  2. * (mindiff > -1) * (maxdiff != 0);
            Fulladjtable(cr, cc) = 1. * (maxdiff == 1);
        }
    }
    arma::uvec finalcols = find(model == 1);
    arma::mat Atable = FullAtable.cols(finalcols);
    arma::mat LBtable = FullLBtable.cols(finalcols);
    arma::mat DtoQtable = FullDtoQtable.cols(finalcols);
    arma::mat adjtable = Fulladjtable.submat(finalcols, finalcols);
    // return FullAtable.cols(finalcols);
    return Rcpp::List::create(
        Rcpp::Named("Atable", Atable), Rcpp::Named("LBtable", LBtable),
        Rcpp::Named("finalcols", finalcols),
        Rcpp::Named("DtoQtable", DtoQtable), Rcpp::Named("adjtable", adjtable));
}

// [[Rcpp::export]]
arma::mat QfromD(unsigned int J, unsigned int K, const arma::mat &DELTA,
                 const arma::mat &DtoQtable)
{
    arma::mat Q(J, K);
    for (unsigned int k = 0; k < K; k++) {
        arma::uvec activek = find(DtoQtable.row(k) == 1);
        arma::mat DELTAk = DELTA.cols(activek);
        for (unsigned int j = 0; j < J; j++) {
            Q(j, k) = arma::max(DELTAk.row(j));
        }
    }
    return Q;
}

// [[Rcpp::export]]
arma::vec permuteAtableIndices(unsigned int nClass, unsigned int K,
                               unsigned int M, unsigned int order,
                               const arma::vec &vv, const arma::vec &perm)
{
    arma::vec vvperm(K);
    arma::vec fullorigperm = arma::linspace(0, nClass - 1, nClass);
    for (unsigned int k = 0; k < K; k++) {
        vvperm(k) = vv(perm(k));
    }
    arma::vec model(nClass);
    arma::vec fullpermindices(nClass);
    for (unsigned int cr = 0; cr < nClass; cr++) {
        arma::vec alpha_r = inv_gen_bijectionvector(K, M, cr);
        double nof0s = 0.;
        for (unsigned int k = 0; k < K; k++) {
            nof0s += 1. * (alpha_r(k) == 0);
        }
        model(cr) = 1. * (nof0s > double(K - order) - 1);
        arma::vec alpha_perm(K);
        fullpermindices(cr) = arma::accu(alpha_r % vvperm);
    }
    arma::uvec finalcols = find(model == 1);
    arma::vec origperm = fullorigperm(finalcols);
    arma::vec reducedpermindices = fullpermindices(finalcols);
    arma::vec permindices(origperm.n_elem);
    for (unsigned int p = 0; p < origperm.n_elem; p++) {
        double origval = origperm(p);
        for (unsigned int pp = 0; pp < origperm.n_elem; pp++) {
            if (origval == reducedpermindices(pp)) {
                permindices(p) = pp;
            }
        }
    }
    return permindices;
}

// [[Rcpp::export]]
double compute_srmr(const arma::rowvec &obs_mean, const arma::mat &obs_cov,
                    const arma::rowvec &est_mean, const arma::mat &est_cov)
{
    double S = 0;
    unsigned int J = obs_mean.n_elem;

    for (unsigned int j1 = 0; j1 < J; j1++) {
        double tmpmeandiff = obs_mean(j1) / sqrt(obs_cov(j1, j1)) -
                             est_mean(j1) / sqrt(est_cov(j1, j1));
        double tmpvardiff =
            (obs_cov(j1, j1) - est_cov(j1, j1)) / obs_cov(j1, j1);
        S += tmpmeandiff * tmpmeandiff + tmpvardiff * tmpvardiff;
    }
    for (unsigned int j1 = 0; j1 < J - 1; j1++) {
        for (unsigned int j2 = j1 + 1; j2 < J; j2++) {
            double tmpcovdiff =
                obs_cov(j1, j2) / sqrt(obs_cov(j1, j1) * obs_cov(j2, j2)) -
                est_cov(j1, j2) / sqrt(est_cov(j1, j1) * est_cov(j2, j2));
            S += tmpcovdiff * tmpcovdiff;
        }
    }
    return sqrt(S / double((J + 1) * J / 2. + J));
}

// [[Rcpp::export]]
Rcpp::List ohoegdm_cpp(const arma::mat &Y, unsigned int K, unsigned int M,
                       unsigned int order, const arma::vec &l0,
                       const arma::vec &l1, 
                       double m0, double bq,
                       double sdMH,
                       unsigned int burnin,
                       unsigned int chain_length = 10000)
{

    // need to allow Ms to vary by j

    unsigned int Malpha = 2;
    unsigned int N = Y.n_rows;
    unsigned int J = Y.n_cols;
    unsigned int nClass = pow(Malpha, K);
    
    unsigned int chain_length_plus_burnin = chain_length + burnin;
    // unsigned int chain_m_burn = chain_length - burnin;
    unsigned int tmburn;
    arma::vec vv = gen_bijectionvector(K, Malpha);

    Rcpp::List designmats = GenerateAtable(nClass, K, Malpha, order);
    arma::mat Atable = Rcpp::as<arma::mat>(designmats[0]);
    arma::mat LBtable = Rcpp::as<arma::mat>(designmats[1]);
    arma::uvec Bindices = Rcpp::as<arma::uvec>(designmats[2]);
    arma::mat qtable = Rcpp::as<arma::mat>(designmats[3]);

    unsigned int P = Atable.n_cols;
    
    // Note: Here we use the full order Atable and then subset out rows needed.
    arma::mat Q_prime = Q_prime_matrix(K, Atable, vv);
    // Q_prime = Q_prime.rows(Bindices); // enable for lower order
    
    // Descales nicely, we hope.
    arma::mat ETA_prime = eta_dina_matrix(Q_prime);
    arma::mat CLtable = CL_gen_invbijection_table(K, Malpha, nClass);
    arma::mat Q = random_Q(J, K);

    arma::mat Q_item_encoded(J, chain_length);
    
    // need to initialize parameters
    arma::vec Tau = arma::zeros<arma::vec>(K);
    arma::vec lambda = arma::ones<arma::vec>(K) * .5;
    arma::vec theta = arma::randn<arma::vec>(N);
    arma::mat KAPPA = kappa_initialize(M, J);
    
    // 
    arma::vec CLASS =
        arma::randi<arma::vec>(N, arma::distr_param(0, nClass - 1));
    
    // arma::vec zerotoMm2=arma::linspace(0,M-2,M-1);
    // for(unsigned int j=0;j<J;j++){
    //   KAPPA.row(j)=2.*zerotoMm2.t();
    // }

    // // arma::mat KAPPA(J,M-1);
    // double intervalbound=3;
    // for(unsigned int j=0;j<J;j++){
    //   for(unsigned int m=0;m<M-1;m++){
    //     KAPPA(j,m)=-intervalbound+(m+1.)*2.*intervalbound/double(M);
    //   }
    // }
    //
    // Rcpp::Rcout << KAPPA  << std::endl;
    // arma::mat DELTA = arma::randi<arma::mat>(J,P,arma::distr_param(0,1));
    // DELTA.col(0)=arma::ones<arma::vec>(J);
    
    // Relate the randomly sampled Q to DELTA through a realization of the
    // the Q prime/ideal matrix under order M.
    
    arma::mat DELTA = q_to_delta(Q, Q_prime, Malpha); 
    
    arma::mat BETA = arma::randu<arma::mat>(J, P);
    double Bplb;
    for (unsigned int j = 0; j < J; j++) {
        arma::rowvec betaj = BETA.row(j);
        arma::rowvec deltaj = DELTA.row(j);
        // arma::rowvec deltaj=DELTA.row(j);
        for (unsigned int p = 0; p < P; p++) {
            double djp = deltaj(p);
            double vforb = djp * l1(p) + (1. - djp) * l0(p);
            double vb = 1. / vforb;
            double mub = 0.;
            mub += 1. * (p == 0) * vb * vforb * m0;
            double sdb = sqrt(vb);
            if (p > 0) {
                // find Bp lowerbound
                // double Bp=betaj(p);
                Bplb =
                    0; // computeLB(nClass,LBtable,Atable,p,betaj,Bp,Bindices);
                if ((Bplb - mub) / sdb > 4.264891) {
                    betaj(p) = Bplb;
                } else {
                    betaj(p) = rTruncNorm_lb(mub, sdb, Bplb);
                }
            } else {
                betaj(p) = R::rnorm(mub, sqrt(vb));
            }
        }
        BETA.row(j) = betaj;
    }
    // double a0(1.),b0(1.);
    
    // arma::vec d0 = dpi0 * arma::ones<arma::vec>(nClass);
    // arma::vec pis = rDirichlet(d0);
    
    
    // Defined between 0 and 1.
    double slipping = R::runif(0, 1);
    double guessing = R::runif(0, 1 - slipping);
    
    double omega(.5);
    
    // must update compute PYa
    Rcpp::List inputs = computePYa(J, M, nClass, Atable, BETA, KAPPA);
    arma::mat ABETA = Rcpp::as<arma::mat>(inputs[0]);
    arma::vec ABETA_sqnorm = Rcpp::as<arma::vec>(inputs[1]);
    arma::cube PY_a = Rcpp::as<arma::cube>(inputs[2]);

    // initialize augmented data from model
    arma::mat Yast(N, J);
    for (unsigned int i = 0; i < N; i++) {
        for (unsigned int j = 0; j < J; j++) {
            double class_i = CLASS(i);
            double Yij = Y(i, j);
            double aiBj = ABETA(j, class_i);
            arma::vec pYij = PY_a.tube(j, class_i);
            Yast(i, j) = rTruncNorm(aiBj, 1., Yij, pYij);
        }
    }

    // ppps output
    /*arma::mat Ymatppp=arma::zeros<arma::mat>(N,J);
    // arma::mat meandiffmat=arma::zeros<arma::mat>(N,J);
    //arma::mat Ytot(N,chain_m_burn);
    //computing observed means and covariance
    arma::rowvec itemmeans = arma::mean(Y);
    arma::mat itemcov = Y.t()*Y;//arma::cov(Y);
    //creating objects for computing ppps on item means and covariances
    arma::rowvec itemmeanppp=arma::zeros<arma::rowvec>(J);
    arma::mat itemcovppp = arma::zeros<arma::mat>(J,J);
    arma::vec srmr(chain_m_burn);*/
    double loglike = 0.;
    
    // Saving output
    arma::cube BETAS(J, P, chain_length);
    arma::mat CLs(N, chain_length);
    arma::mat PIs(nClass, chain_length);
    arma::vec omegas(chain_length);
    // arma::vec omegaTaus(chain_length);
    arma::cube KAPPAs = arma::cube(J, M - 1, chain_length);
    // arma::mat mKAPPAs=arma::cube(J,M-1);
    arma::vec mtheta = arma::zeros<arma::vec>(N);
    // arma::mat D_tab=arma::zeros<arma::mat>(J,P);
    arma::mat Taus(K, chain_length);
    arma::mat lambdas(K, chain_length);
    arma::mat Q_tab = arma::zeros<arma::mat>(J, K);
    
    arma::vec m2lls(chain_length);
    arma::vec SLIP(chain_length);
    arma::vec GUESS(chain_length);
    
    arma::mat DELTA_tab = arma::zeros<arma::mat>(J, P);
    arma::vec MHsum = arma::zeros<arma::vec>(J);
    arma::vec MHaccept = arma::zeros<arma::vec>(J);
    
    // Start Markov chain
    for (unsigned int t = 0; t < chain_length_plus_burnin; t++) {
        omega = parm_update_nomiss(
            N, J, K, nClass, M, Y, Yast, BETA, KAPPA, CLASS, theta, lambda, Tau,
            Q, 
            DELTA,
            Q_prime,
            ETA_prime,
            slipping,
            guessing,
            omega, vv, CLtable, Atable, LBtable, Bindices, qtable, P, l1, m0,
            l0,
            bq, ABETA, ABETA_sqnorm, PY_a, MHaccept, sdMH, loglike);

        
        // update inputs
        if (t > burnin - 1) {
            tmburn = t - burnin;
            /*
          arma::mat Ysim=sim_slcm(N,J,M,nClass,CLASS,Atable,BETA,KAPPA);
          arma::rowvec sim_item_means =arma::mean(Ysim);
          arma::mat sim_item_cov = Ysim.t()*Ysim;//arma::cov(Ysim);
          for(unsigned int j1=0;j1<J;j1++){
            itemmeanppp(j1)+= 1.*(itemmeans(j1)>sim_item_means(j1));
            itemcovppp(j1,j1) += 1.*(itemcov(j1,j1)>sim_item_cov(j1,j1));
          }
          for(unsigned int j1=0;j1<J-1;j1++){
            for(unsigned int j2=j1+1;j2<J;j2++){
              itemcovppp(j1,j2) += 1.*(itemcov(j1,j2)>sim_item_cov(j1,j2));
            }
          }
          srmr(tmburn) =
          compute_srmr(itemmeans,itemcov,sim_item_means,sim_item_cov);

          // Ytot.col(tmburn)=arma::sum(Ysim,1);

          for(unsigned int i=0;i<N;i++){
            for(unsigned int j=0;j<J;j++){
              Ymatppp(i,j)+=1.*(Ysim(i,j)==Y(i,j));
              // meandiffmat(i,j)+=Y(i,j)-Ysim(i,j);
            }
          }*/
            m2lls(tmburn) = -2. * loglike; // slcm_m2LL(N,J,M,nClass,Y,pis,PY_a);
            Q_item_encoded.col(tmburn) = Q * vv; 
            SLIP(tmburn) = slipping;
            GUESS(tmburn) = guessing;
            BETAS.slice(tmburn) = BETA;
            KAPPAs.slice(tmburn) = KAPPA;
            CLs.col(tmburn) = CLASS;
            // mKAPPAs   = KAPPA/double(tmburn+1.)
            // +double(tmburn/(tmburn+1.))*mKAPPAs;
            mtheta = theta / double(tmburn + 1.) +
                     double(tmburn / (tmburn + 1.)) * mtheta;
            // D_tab               +=DELTA;
            Q_tab += Q;
            DELTA_tab += DELTA;
            
            omegas(tmburn) = omega;
            MHsum += MHaccept;
            // MHaccept=arma::zeros<arma::vec>(J);
            Taus.col(tmburn) = Tau;
            lambdas.col(tmburn) = lambda;
        }
    }
    
    Rcpp::List estimates = 
        Rcpp::List::create(
            Rcpp::Named("mtheta", mtheta),
            Rcpp::Named("QS", Q_tab / chain_length),
            Rcpp::Named("deltas", DELTA_tab / chain_length)
        );

    Rcpp::List chain = 
        Rcpp::List::create(
        Rcpp::Named("betas", BETAS),
        Rcpp::Named("guessing", GUESS),
        Rcpp::Named("slipping", SLIP),
        Rcpp::Named("kappas", KAPPAs), // CHANGE IN THE CODE ABOVE
        Rcpp::Named("taus", Taus),  // CHANGE IN THE CODE ABOVE
        Rcpp::Named("m2lls", m2lls),
        Rcpp::Named("omegas", omegas),
        Rcpp::Named("lambdas", lambdas),
        Rcpp::Named("classes", CLs)
        );

    Rcpp::List recovery = 
      Rcpp::List::create(
        Rcpp::Named("Q_item_encoded", Q_item_encoded),
        Rcpp::Named("MHsum", MHsum / chain_length)
      );
    
    return Rcpp::List::create(
        Rcpp::Named("estimates", estimates),
        Rcpp::Named("chain", chain),
        Rcpp::Named("recovery", recovery)
    );
}

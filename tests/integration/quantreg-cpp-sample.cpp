#include <quantreg-cpp>

int main() {
  try {
    // engel dataset
    arma::mat data;
    data.load("./data/engel.csv", arma::csv_ascii);
    data.shed_cols(0, 0);
    data.shed_rows(0, 0);
    arma::vec y = data.col(1);
    arma::mat X = data(arma::span::all, arma::span(0, 0));
    double tau = 0.5;

    arma::vec coefficients;
    arma::vec residuals;
    std::tie(coefficients, residuals) = quantreg::qr(X, y, tau);

    coefficients.t().print(std::cout, "coef=");
    // std::cout.precision(10);
    // coefficients.t().raw_print("coef=");

  } catch (const std::runtime_error& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Error: unknown exception" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

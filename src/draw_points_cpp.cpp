#include <Rcpp.h>
#include <Rinternals.h>

//' Draw points.
//'
//' @param total_points Total number of points to draw.
//' @param draw_points_lambda Function to draw the points.
//' @param lambda_max Upper bound to lambda.
//'
//' @useDynLib compp
//' @import Rcpp
// [[Rcpp::export]]
Rcpp::NumericMatrix draw_points_cpp(int total_points,
                                    Rcpp::Function draw_points_lambda,
                                    double lambda_max) {
  Rcpp::NumericMatrix points(Rcpp::no_init(total_points, 2));

  Rcpp::NumericMatrix lambda_draws(draw_points_lambda());
  int index(0);

  const int simultaneous_draws(lambda_draws.nrow());

  for(int i(0); i < total_points; ++i) {
    while(lambda_draws(index, 2) < unif_rand() * lambda_max) {
      if(index >= simultaneous_draws - 1) {
        lambda_draws = draw_points_lambda();
        index = 0;
      } else {
        ++index;
      }
    }
    points(i, 0) = lambda_draws(index, 0);
    points(i, 1) = lambda_draws(index, 1);
    if(index >= simultaneous_draws - 1) {
      lambda_draws = draw_points_lambda();
      index = 0;
    } else {
      ++index;
    }
  }
  return points;
}

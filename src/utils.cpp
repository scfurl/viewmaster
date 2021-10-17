// we only include RcppArrayFire.h which pulls Rcpp.h in for us
#include "RcppArrayFire.h"
#include <arrayfire.h>
#include <math.h>
#include <stdio.h>
#include <af/util.h>
#include <string>
#include <vector>
#include "mnist_common.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppFire so that the build process will know what to do
//
// [[Rcpp::depends(RcppArrayFire)]]

// RcppArrayFire needs C++11
// add the following comment when you export your
// C++ function to R via Rcpp::SourceCpp()
// [[Rcpp::plugins(cpp11)]]

// simple example of creating two matrices and
// returning the result of an operation on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//

using namespace af;

Rcpp::List return_mnist(int perc, bool verbose = false) {
  array train_images, train_labels;
  array test_images, test_labels;
  int num_train, num_test, num_classes;
  // Load mnist data
  float frac = (float)(perc) / 100.0;
  setup_mnist<false>(&num_classes, &num_train, &num_test, train_images,
                     test_images, train_labels, test_labels, frac);
  int feature_length = train_images.elements() / num_train;
  array train_feats  = moddims(train_images, feature_length, num_train);
  array test_feats   = moddims(test_images, feature_length, num_test);
  // Get training parameters
  array mu, sig2;
  if(verbose){
    std::cerr << "Train feature dims:" << std::endl;
    std::cerr << train_feats.dims() << std::endl;
    std::cerr << "Test feature dims:" << std::endl;
    std::cerr << test_feats.dims() << std::endl;
    std::cerr << "Train labels dims:" << std::endl;
    std::cerr << train_labels.dims() << std::endl;
    std::cerr << "Test labels dims:" << std::endl;
    std::cerr << test_labels.dims() << std::endl;
    std::cerr << "Num classes:" << std::endl;
    std::cerr << num_classes << std::endl;
  }
  Rcpp::List rl = Rcpp::List::create(train_feats, test_feats, train_labels, test_labels);
  return(rl);
}

//' @export
// [[Rcpp::export]]
Rcpp::List get_mnist(int perc = 80, bool verbose = true) {
  // int device   = argc > 1 ? atoi(argv[1]) : 0;
  // bool console = argc > 2 ? argv[2][0] == '-' : false;
  // int perc     = argc > 3 ? atoi(argv[3]) : 60;
  af::setDevice(0);
  std::string info_string = af::infoString();
  std::cerr << info_string;
  return(return_mnist(perc, verbose));
  // try {
  //     af::setDevice(0);
  //     af::info();
  //     naive_bayes_demo(perc);
  // } catch (af::exception &ae) { std::cerr << ae.what() << std::endl; }
}
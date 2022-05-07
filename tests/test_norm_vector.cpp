//
// Created by evgen on 19.02.2022.
//

#include "gtest/gtest.h"
#include "my_project/vectors/Norm.hpp"
#include <vector>

class VECTOR_NORM_TESTS: public ::testing::Test {
protected:
    void SetUp() override { };
    std::vector<double> v1 = {1, 1, 1, 1, 1, 1};
    std::vector<double> v2 = {-1, -2, -3, -4, -5, -6, -7, 10, -11};
    std::vector<double> v3 = {1.1, 0, 0, 0, 0, 0};
    std::vector<double> v4 = {-1.1, 0, 0, 0, 0, 0};
    std::vector<double> v5 = {-1.1, -1.1, 2.2, 2.3, -2.3, -2.2};
    std::vector<double> v6 = {0.};
    double tolerance = 0.00001;
};

TEST_F(VECTOR_NORM_TESTS, FIRST_NORM){
    ASSERT_NEAR(norm(v1, NormType::FirstNorm), 6, tolerance) << "First norm dont work properly";
    ASSERT_NEAR(norm(v2, NormType::FirstNorm), 49, tolerance) << "First norm dont work properly";
    ASSERT_NEAR(norm(v3, NormType::FirstNorm), 1.1, tolerance) << "First norm dont work properly";
    ASSERT_NEAR(norm(v4, NormType::FirstNorm), 1.1, tolerance) << "First norm dont work properly";
    ASSERT_NEAR(norm(v5, NormType::FirstNorm), 11.2, tolerance) << "First norm dont work properly";
    ASSERT_NEAR(norm(v6, NormType::FirstNorm), 0, tolerance) << "First norm dont work properly";
}

TEST_F(VECTOR_NORM_TESTS, FIRST_NORM_NEW){
    ASSERT_NEAR(Norm<NormType::FirstNorm>(v1), 6, tolerance) << "First norm dont work properly";
    ASSERT_NEAR(Norm<NormType::FirstNorm>(v2), 49, tolerance) << "First norm dont work properly";
    ASSERT_NEAR(Norm<NormType::FirstNorm>(v3), 1.1, tolerance) << "First norm dont work properly";
    ASSERT_NEAR(Norm<NormType::FirstNorm>(v4), 1.1, tolerance) << "First norm dont work properly";
    ASSERT_NEAR(Norm<NormType::FirstNorm>(v5), 11.2, tolerance) << "First norm dont work properly";
    ASSERT_NEAR(Norm<NormType::FirstNorm>(v6), 0, tolerance) << "First norm dont work properly";
}

TEST_F(VECTOR_NORM_TESTS, VECTOR_SECOND_NORM){
    ASSERT_NEAR(norm(v1, NormType::SecondNorm), std::sqrt(6), tolerance) << "Second norm dont work properly";
    ASSERT_NEAR(norm(v2, NormType::SecondNorm), std::sqrt(361), tolerance) << "Second norm dont work properly";
    ASSERT_NEAR(norm(v3, NormType::SecondNorm), 1.1, tolerance) << "Second norm dont work properly";
    ASSERT_NEAR(norm(v4, NormType::SecondNorm), 1.1, tolerance) << "Second norm dont work properly";
    ASSERT_NEAR(norm(v5, NormType::SecondNorm), std::sqrt(22.68), tolerance) << "Second norm dont work properly";
    ASSERT_NEAR(norm(v6, NormType::SecondNorm), std::sqrt(0), tolerance) << "Second norm dont work properly";
};

TEST_F(VECTOR_NORM_TESTS, VECTOR_SECOND_NORM_NEW){
    ASSERT_NEAR(Norm<NormType::SecondNorm>(v1), std::sqrt(6), tolerance) << "Second norm dont work properly";
    ASSERT_NEAR(Norm<NormType::SecondNorm>(v2), std::sqrt(361), tolerance) << "Second norm dont work properly";
    ASSERT_NEAR(Norm<NormType::SecondNorm>(v3), 1.1, tolerance) << "Second norm dont work properly";
    ASSERT_NEAR(Norm<NormType::SecondNorm>(v4), 1.1, tolerance) << "Second norm dont work properly";
    ASSERT_NEAR(Norm<NormType::SecondNorm>(v5), std::sqrt(22.68), tolerance) << "Second norm dont work properly";
    ASSERT_NEAR(Norm<NormType::SecondNorm>(v6), std::sqrt(0), tolerance) << "Second norm dont work properly";
};

TEST_F(VECTOR_NORM_TESTS, VECTOR_INF_NORM){
    ASSERT_NEAR(norm(v1, NormType::InfNorm), 1, tolerance) <<  "Inf norm dont work properly";
    ASSERT_NEAR(norm(v2, NormType::InfNorm), 11, tolerance) <<  "Inf norm dont work properly";
    ASSERT_NEAR(norm(v3, NormType::InfNorm), 1.1, tolerance) <<  "Inf norm dont work properly";
    ASSERT_NEAR(norm(v4, NormType::InfNorm), 1.1, tolerance) <<  "Inf norm dont work properly";
    ASSERT_NEAR(norm(v5, NormType::InfNorm), 2.3, tolerance) <<  "Inf norm dont work properly";
    ASSERT_NEAR(norm(v6, NormType::InfNorm), 0, tolerance) <<  "Inf norm dont work properly";
};

TEST_F(VECTOR_NORM_TESTS, VECTOR_INF_NORM_NEW){
    ASSERT_NEAR(Norm<NormType::InfNorm>(v1), 1, tolerance) <<  "Inf norm dont work properly";
    ASSERT_NEAR(Norm<NormType::InfNorm>(v2), 11, tolerance) <<  "Inf norm dont work properly";
    ASSERT_NEAR(Norm<NormType::InfNorm>(v3), 1.1, tolerance) <<  "Inf norm dont work properly";
    ASSERT_NEAR(Norm<NormType::InfNorm>(v4), 1.1, tolerance) <<  "Inf norm dont work properly";
    ASSERT_NEAR(Norm<NormType::InfNorm>(v5), 2.3, tolerance) <<  "Inf norm dont work properly";
    ASSERT_NEAR(Norm<NormType::InfNorm>(v6), 0, tolerance) <<  "Inf norm dont work properly";
};

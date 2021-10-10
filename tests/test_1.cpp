#include <catch2/catch.hpp>

#include <nasoq/nasoq.h>
#include <nasoq/nasoq_eigen.h>

#include <fstream>
#include <nlohmann/json.hpp>

#include <iostream>

using Json = nlohmann::json;


struct TestData
{
    Eigen::SparseMatrix<double> H;
    Eigen::VectorXd q;
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b;
    Eigen::SparseMatrix<double> C;
    Eigen::VectorXd d;
    Eigen::VectorXd x;
    Eigen::VectorXd y;
    Eigen::VectorXd z;
};


void load_json(Eigen::VectorXd& out_v, const Json& json) {
    int num_elements = static_cast<int>(json.size());
    out_v.resize(num_elements);
    for(int i = 0; i < num_elements; ++i) {
        out_v[i] = json[i];
    }
}


void load_json(Eigen::SparseMatrix<double>& out_M, const Json& json) {
    int rows = json["rows"].get<int>();
    int cols = json["cols"].get<int>();
    bool is_compressed = !json.contains("inner_nnzs");
    REQUIRE(is_compressed);

    out_M.resize(rows, cols);
    int outer_size = out_M.outerSize();

    const Json& json_outer_starts = json["outer_starts"];
    for(int outer_ind = 0; outer_ind <= outer_size; ++outer_ind) {
        out_M.outerIndexPtr()[outer_ind] = json_outer_starts[outer_ind].get<int>();
    }

    int values_size = out_M.outerIndexPtr()[outer_size];
    out_M.reserve(values_size);
    out_M.resizeNonZeros(values_size);

    const Json& json_inner_inds = json["inner_indices"];
    const Json& json_values = json["values"];
    for(int i = 0; i < values_size; ++i) {
        out_M.innerIndexPtr()[i] = json_inner_inds[i].get<int>();
        out_M.valuePtr()[i] = json_values[i].get<double>();
    }
}


void load_json(TestData& out_data, const Json& json) {
    load_json(out_data.H, json["H"]);
    load_json(out_data.q, json["q"]);
    load_json(out_data.A, json["A"]);
    load_json(out_data.b, json["b"]);
    load_json(out_data.C, json["C"]);
    load_json(out_data.d, json["d"]);
    load_json(out_data.x, json["x"]);
    load_json(out_data.y, json["y"]);
    load_json(out_data.z, json["z"]);
}


TEST_CASE("nasoq::quadprog test 1") {
    constexpr int n = 41;

    std::ifstream i(TEST_DIR "/1.json");
    Json json;
    i >> json;
    TestData test_data;
    load_json(test_data, json);

    auto* qs = new nasoq::QPSettings();

    Eigen::VectorXd x;
    Eigen::VectorXd y;
    Eigen::VectorXd z;
    int ret = nasoq::quadprog(test_data.H, test_data.q, test_data.A, test_data.b, test_data.C, test_data.d, x, y, z, qs);
    REQUIRE(ret > 0);

    REQUIRE(x.isApprox(test_data.x));
    REQUIRE(y.isApprox(test_data.y));
    REQUIRE(z.isApprox(test_data.z));

}

#pragma once
#include <memory>
#include <vector>
#include <string>


namespace utils {
    constexpr double EPS = 1e-16;

    inline double e(int i, int j) {
        return (i == j);
    }

    inline double U(int i, int j) {
        if (i > j) return 0;
        else if (i == j) return 1;
        else return 1. / (i + j + 1);
    }

    inline double hilbert(int i, int j) {
        return 1. / (i + j + 1);
    }

    inline double absDiff(int i, int j) {
        return abs(i - j);
    }

    inline double revAbsDiff(int i, int j) {
        return 1. / (absDiff(i, j) + 1);
    }

    inline bool isPositive(double value) {
        return (value > EPS);
    }

    inline bool isEqual(double lhs, double rhs) {
        return (abs(lhs - rhs) <= EPS);
    }

    std::vector<std::string> split(const std::string& text, char sep) {
        std::vector<std::string> tokens;
        std::size_t start = 0, end = 0;
        while ((end = text.find(sep, start)) != std::string::npos) {
            tokens.push_back(text.substr(start, end - start));
            start = end + 1;
        }
        tokens.push_back(text.substr(start));
        return tokens;
    };
};
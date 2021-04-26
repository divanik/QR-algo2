#pragma once

#include<set>
#include<utility>
#include<vector>

namespace QR_algorithm {


template <typename T>
class Shift_splitter {

public:
    using Iterator = std::set<std::pair<int, int>>::iterator;
    using Const_iterator = std::set<std::pair<int, int>>::const_iterator;

    Shift_splitter(size_t lef_, size_t rig_);

    void fill_splitter(int ind);

    void split_segs(size_t lef, size_t rig);

    void flush_buffer(bool);

    Iterator begin();

    Iterator end();

    Const_iterator begin() const;

    Const_iterator end() const;

private:
    std::vector<int> splitter;
    std::set<std::pair<int, int>> segs;
    std::set<std::pair<int, int>> to_ins;
    std::set<std::pair<int, int>> to_del;
    MatrixX<T>* matr;
};

template<typename T>
Shift_splitter<T>::Shift_splitter(size_t lef_, size_t rig_) {
    segs.insert({lef_, rig_});
}
template<typename T>
void Shift_splitter<T>::fill_splitter(int ind) {
    splitter.push_back(ind);
}

template<typename T>
void Shift_splitter<T>::split_segs(size_t lef, size_t rig) {
    if (splitter.size()) {
        to_del.insert({lef, rig});
        int cur_lef = lef;
        for (auto& x : splitter) {
            if (cur_lef != x) {
                to_ins.insert({cur_lef, x});
            }
            cur_lef = x + 1;
        }
        if (cur_lef != rig) {
            to_ins.insert({cur_lef, rig});
        }
        splitter.clear();
    }
}

template <typename T>
void Shift_splitter<T>::flush_buffer(bool pseudo_shur = false) {
    for (auto& x : to_del) {
        segs.erase(x);
    }
    to_del.clear();

    for (auto& x : to_ins) {
        if (pseudo_shur) {
            if (x.second == x.first + 1) {
                T screw_trace = matr(x.first, x.first) - matr(x.second, x.second);
                T det = screw_trace * screw_trace + 4 * matr(x.first, x.second) * matr(x.second, x.first);
                if (det >= T(0)) {
                    return det;
                }
            } else {
                segs.insert(x);
            }
        }
    }
    to_ins.clear();
}

template <typename T>
Shift_splitter<T>::Iterator Shift_splitter<T>::begin() {
    return segs.begin();
}

template <typename T>
Shift_splitter<T>::Iterator Shift_splitter<T>::end() {
    return segs.end();
}

template <typename T>
Shift_splitter<T>::Const_iterator Shift_splitter<T>::begin() const {
    return segs.begin();
}

template <typename T>
Shift_splitter<T>::Const_iterator Shift_splitter<T>::end() const{
    return segs.end();
}



}
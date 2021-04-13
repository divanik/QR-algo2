#pragma once

#include<set>
#include<utility>
#include<vector>

namespace QR_algorithm {

class Shift_splitter {

public:
    using Iterator = std::set<std::pair<int, int>>::iterator;
    using Const_iterator = std::set<std::pair<int, int>>::const_iterator;

    Shift_splitter(size_t lef_, size_t rig_);

    void fill_splitter(int ind);

    void split_segs(size_t lef, size_t rig);

    void flush_buffer();

    Iterator begin();

    Iterator end();

    Const_iterator begin() const;

    Const_iterator end() const;

private:
    std::vector<int> splitter;
    std::set<std::pair<int, int>> segs;
    std::set<std::pair<int, int>> to_ins;
    std::set<std::pair<int, int>> to_del;
};

}
#include "shift_splitter.h"

#include<set>
#include<utility>
#include<vector>

namespace QR_algorithm {

    Shift_splitter::Shift_splitter(size_t lef_, size_t rig_) {
        segs.insert({lef_, rig_});
    }

    void Shift_splitter::fill_splitter(int ind) {
        splitter.push_back(ind);
    }

    void Shift_splitter::split_segs(size_t lef, size_t rig) {
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

    void Shift_splitter::flush_buffer() {
        for (auto& x : to_del) {
            segs.erase(x);
        }
        to_del.clear();

        for (auto& x : to_ins) {
            segs.insert(x);
        }
        to_ins.clear();
    }

    Shift_splitter::Iterator Shift_splitter::begin() {
        return segs.begin();
    }

    Shift_splitter::Iterator Shift_splitter::end() {
        return segs.end();
    }

    Shift_splitter::Const_iterator Shift_splitter::begin() const {
        return segs.begin();
    }

    Shift_splitter::Const_iterator Shift_splitter::end() const{
        return segs.end();
    }

}
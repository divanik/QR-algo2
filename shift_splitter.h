#include<set>
#include<utility>
#include<vector>

class Shift_splitter {

public:
    Shift_splitter(size_t lef_, size_t rig_) {
        segs.insert({lef_, rig_});
    }

    void fill_splitter(int ind) {
        splitter.push_back(ind);
    }

    void split_segs(size_t lef, size_t rig) {
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

    void flush_buffer() {
        for (auto& x : to_del) {
            segs.erase(x);
        }
        to_del.clear();

        for (auto& x : to_ins) {
            segs.insert(x);
        }
        to_ins.clear();
    }

    std::set<std::pair<int, int>>::iterator begin() {
        return segs.begin();
    }

    std::set<std::pair<int, int>>::iterator end() {
        return segs.end();
    }

private:
    std::vector<int> splitter;
    std::set<std::pair<int, int>> segs;
    std::set<std::pair<int, int>> to_ins;
    std::set<std::pair<int, int>> to_del;
};
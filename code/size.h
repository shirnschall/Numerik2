//
// Created by Sebastian Hirnschall on 18.01.20.
//

#ifndef AUFGABE1_SIZE_H
#define AUFGABE1_SIZE_H

namespace linag {
    class Size {
    public:
        int rows, cols;



    };

    bool operator==(const Size& lhs, const Size& rhs)
    {
        return lhs.rows == rhs.rows && lhs.cols == rhs.cols;
    }

    bool operator!=(const Size& lhs, const Size& rhs)
    {
        return !(lhs == rhs);
    }
}

#endif //AUFGABE1_SIZE_H

//
// Created by Sebastian Hirnschall on 17.01.20.
//

#ifndef AUFGABE1_VECTOR_H
#define AUFGABE1_VECTOR_H

namespace linag {

    template<typename T>
    class Vector {
    private:
        int length;
        T *data;
    public:
        Vector();
        ~Vector();

        Vector(Vector<T>);
        Vector<T> &operator=(const Vector<T> &);
        Vector(int length);

        T& at(int);
        const T &at(int) const;

        int dim();
    };

    template<typename T>
    std::ostream& operator<<(std::ostream& output,const Vector<T>& x);


}


#endif //AUFGABE1_VECTOR_H

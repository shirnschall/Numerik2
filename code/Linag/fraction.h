#ifndef CODE_FRACTION_H
#define CODE_FRACTION_H


#include <iostream>
#include <cassert>
#include <cmath>

using std::cout;
namespace linag {
    class Fraction {
    private:
        int p;  //zaehler
        int q;  //nenner
    public:
        Fraction(int p = 0, int q = 1);  //standardkonstruktor mit default werten
        Fraction(const Fraction &);  //kopierkonstruktor
        Fraction(double a); //typekast von double
        ~Fraction();    //destruktor

        Fraction &operator=(const Fraction &);

        const Fraction operator-() const;
        Fraction& operator-=(const Fraction& rhs);
        Fraction& operator+=(const Fraction& rhs);
        Fraction& operator*=(const Fraction& rhs);
        Fraction& operator/=(const Fraction& rhs);

        operator double() const;  //typecast auf double
        operator int() const;

        void reduce();  //kuerzen
        int getGGT(int a, int b);    //euklid

        int getNumerator() const { return p; };

        int getDenominator() const { return q; };

        void setNumerator(int n) { p = n; };

        void setDenominator(int n) {
            assert(n != 0);
            if (n > 0) {
                q = n;
            } else {
                q = -n;
                p = -p;
            }
        }

        void print() {
            cout << p << "/" << q << "\n";
        }
    };

    std::ostream &operator<<(std::ostream &, const Fraction &);

    const Fraction operator+(const Fraction &, const Fraction &);

    const Fraction operator-(const Fraction &, const Fraction &);

    const Fraction operator*(const Fraction &, const Fraction &);

    const Fraction operator/(const Fraction &, const Fraction &);

}


linag::Fraction::Fraction(int p, int q) {
    assert(q!=0);
    if(q < 0)
    {
        this->p = -p;
        this->q = -q;
    }else {
        this->p = p;
        this->q = q;
    }
    reduce();
}

linag::Fraction::Fraction(const Fraction &fraction):
        Fraction(fraction.getNumerator(),fraction.getDenominator()) {}

linag::Fraction::Fraction(double x) {
    //approximation als kettenbruch
    double eps = 1e-14; //fehler
    double a = std::floor(x);
    double y = (double)1/(x-a);
    double p0 = 1;
    double p1 = a;
    double q0 = 0;
    double q1 = 1;
    double tmp;

    do{
        a = std::floor(y);
        y = (double)1/(y-a);
        tmp = p1;
        p1 = p1*a + p0;
        p0 = tmp;

        tmp=q1;
        q1 = q1*a + q0;
        q0 = tmp;
    }while(std::fabs((double)p0/q0 - (double)p1/q1) > eps);

    setNumerator((int)p0);
    setDenominator((int)q0);
    reduce();
}

linag::Fraction::~Fraction() {}

linag::Fraction& linag::Fraction::operator=(const Fraction &fraction) {
    this->p = fraction.getNumerator();
    this->q = fraction.getDenominator();

    return *this;
}

const linag::Fraction linag::Fraction::operator-() const{
    return Fraction(-getNumerator(),getDenominator());
}

linag::Fraction::operator double() const{
    return (double)getNumerator()/getDenominator();
}

linag::Fraction::operator int() const{
    //rest bei division durch 1
    double a = std::fmod((double)getNumerator()/getDenominator(),1);

    if((std::fabs(std::fabs(a) - .5) < 1e-9) || (std::fabs(a) > .5)){
        if(a>0)
            return (int)(double)getNumerator()/getDenominator() + 1;
        else
            return (int)(double)getNumerator()/getDenominator() - 1;
    }else{
        return (int)(double)getNumerator()/getDenominator();
    }
}

void linag::Fraction::reduce() {
    int ggt = getGGT(getNumerator(),getDenominator());
    if(ggt != 1){
        setNumerator(getNumerator()/ggt);
        setDenominator(getDenominator()/ggt);
    }
}

int linag::Fraction::getGGT(int a, int b) {
    //euklid
    a = std::abs(a);
    b = std::abs(b);
    if(a == 0)
        return b;

    while(b != 0){
        if(a>b)
            a -= b;
        else
            b -= a;
    }
    return a;
}

const linag::Fraction linag::operator+(const linag::Fraction& a,const linag::Fraction& b){
    linag::Fraction tmp(a.getNumerator()*b.getDenominator()+b.getNumerator()*a.getDenominator(),
                 a.getDenominator()*b.getDenominator());
    return tmp;
}

const linag::Fraction linag::operator-(const linag::Fraction& a,const linag::Fraction& b){
    return a + (-b);
}

const linag::Fraction linag::operator*(const linag::Fraction& a,const linag::Fraction& b){
    linag::Fraction tmp(a.getNumerator()*b.getNumerator(),
                 a.getDenominator()*b.getDenominator());
    return tmp;
}

const linag::Fraction linag::operator/(const linag::Fraction& a,const linag::Fraction& b){
    assert(std::fabs((double)b.getNumerator()/b.getDenominator())>1e-9);
    linag::Fraction tmp(b.getDenominator(),b.getNumerator());
    return a * tmp;
}

std::ostream& operator<<(std::ostream& output,const linag::Fraction& fraction){
    return output << fraction.getNumerator() << "/" << fraction.getDenominator();
}

linag::Fraction& linag::Fraction::operator-=(const linag::Fraction& rhs){
    double den = getDenominator()*rhs.getDenominator();

    setNumerator(getNumerator()*rhs.getDenominator()-rhs.getNumerator()*getDenominator());
    setDenominator(den);
    return *this;
}

linag::Fraction& linag::Fraction::operator+=(const linag::Fraction& rhs){
    double den = getDenominator()*rhs.getDenominator();

    setNumerator(getNumerator()*rhs.getDenominator()+rhs.getNumerator()*getDenominator());
    setDenominator(den);
    return *this;
}

linag::Fraction& linag::Fraction::operator*=(const linag::Fraction& rhs){
    setDenominator(getDenominator()*rhs.getDenominator());
    setNumerator(getNumerator()*rhs.getNumerator());
    return *this;
}

linag::Fraction& linag::Fraction::operator/=(const linag::Fraction& rhs){
    setDenominator(getDenominator()*rhs.getNumerator());
    setNumerator(getNumerator()*rhs.getDenominator());
    return *this;
}




#endif //CODE_FRACTION_H

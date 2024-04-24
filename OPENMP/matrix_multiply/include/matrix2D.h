class matrix2D
{
private:
    int *mat;
    int n;
    int m;

public:
    matrix2D(int n, int m);
    void set(int x, int y, int value);
    int get(int x, int y);
    int get_n();
    int get_m();
};
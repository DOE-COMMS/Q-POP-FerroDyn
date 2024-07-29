#include <iostream>

using namespace std;

int main()
{
    int nx = 2;
    int ny = 2;
    int nz = 40;

    int n = nx * ny * nz;

    int i, j, k;
    int i1, j1, k1;
    int i2, j2, k2;
    int i3, j3, k3;

    for (int id = 0; id < n; id++)
    {
        i = id / ((ny + 1) * (nz + 1));
		j = (id - i * ((ny + 1) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny + 1) * (nz + 1)) - j * (nz + 1);

        i1 = id / ((ny + 1) * (nz + 1));
        j1 = (id - i1 * ((ny + 1) * (nz + 1))) / (nz + 1);
        k1 = id - i1 * ((ny + 1) * (nz + 1)) - j1 * (nz + 1);

        i2 = id / ((ny) * (nz + 1));
        j2 = (id - i2 * ((ny) * (nz + 1))) / (nz + 1);
        k2 = id - i2 * ((ny) * (nz + 1)) - j2 * (nz + 1);

        i3 = id / ((ny + 1) * (nz));
        j3 = (id - i3 * ((ny + 1) * (nz))) / nz;
        k3 = id - i3 * ((ny + 1) * (nz)) - j3 * nz;

        cout << "(" << i << ", " << j << ", " << k << "); ";
        if (i < nx)
            cout << "i (" << i1 << ", " << j1 << ", " << k1 << "); ";
        if (j < ny)
            cout << "j (" << i2 << ", " << j2 << ", " << k2 << "); ";
        if (k < nz)
            cout << "k (" << i3 << ", " << j3 << ", " << k3 << ")" << endl;
    }
}
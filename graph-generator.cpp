#include <cstdio>
#include <iostream>

using namespace std;

int main()
{
    ios::sync_with_stdio(false);
    int n;
    cin >> n;
    srand(time(NULL));
    cout << n << "\n";
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            cout << rand() % 200 << " ";
    cout << "\n";
    return 0;
}
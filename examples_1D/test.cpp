#include <iostream>
#include <vector>
#include <algorithm> // std::min_element

int main() {
    std::vector<int> numbers = {10, 20, 5, 15, 30};

    // std::min_elementを使って最小の要素を見つける
    auto minElement = std::min_element(numbers.begin(), numbers.end());

    
    int min = *minElement;
    std::cout << min;

    return 0;
}

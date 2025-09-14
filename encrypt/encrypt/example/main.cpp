#include <iostream>
#include "../include/encrypt.hpp"

int main() {
	encryp(0x1);
	std::cout << ("unencrypted string") << std::endl;
	std::cout << encrypt("encrypted string") << std::endl;
	std::cin.get();
	return 0;
}
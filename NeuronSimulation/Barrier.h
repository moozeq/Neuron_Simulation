#pragma once
class Barrier
{
	unsigned channelsIndexFrom;
	unsigned channelsIndexTo;

public:
	Barrier();
	bool checkCollision(float* coordinates);
};


#pragma once
class Util
{
public:
	Util();
	~Util();

	static double TestFunction(double x);

	static bool AtomBalance(double feedStreamMoles[], double prodStreamMoles[], double& CError, double &HError);
	static double brents_fun(std::function<double(double)> f, double lower_bound, double upper_bound, double TOL, int MAX_ITER);

};


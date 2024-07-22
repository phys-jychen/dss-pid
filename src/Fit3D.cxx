#include "Fit3D.h"

vector<TVector3> points;
vector<Double_t> weights;

Fit3D::Fit3D(const vector<Double_t>& hit_x, const vector<Double_t>& hit_y, const vector<Double_t>& hit_z, const vector<Double_t>& hit_energy, const Int_t& nhits)
{
    if (nhits >= 3)
    {
        for (Int_t i = 0; i < nhits; ++i)
        {
            points.emplace_back(hit_x.at(i), hit_y.at(i), hit_z.at(i));
            weights.emplace_back(hit_energy.at(i));
        }

        Math::Minimizer* minimiser = Math::Factory::CreateMinimizer("Minuit2", "Migrad");
        minimiser->SetMaxFunctionCalls(1000000);
        minimiser->SetMaxIterations(100000);
        minimiser->SetTolerance(0.001);

        Math::Functor f(&SumDist2, 6);
        minimiser->SetFunction(f);

        minimiser->SetVariable(0, "x0", 0, 0.1);
        minimiser->SetVariable(1, "y0", 0, 0.1);
        minimiser->SetVariable(2, "z0", 0, 0.1);
        minimiser->SetVariable(3, "dx", 0, 0.001);
        minimiser->SetVariable(4, "dy", 0, 0.001);
        minimiser->SetVariable(5, "dz", 1, 0.001);

        minimiser->Minimize();

        const Double_t* results = minimiser->X();
        TVector3 linePoint(results[0], results[1], results[2]);
        TVector3 lineDir(results[3], results[4], results[5]);
        lineDir = lineDir.Unit();

        direction_x = results[3];
        direction_y = results[4];
        direction_z = results[5];

        Double_t d2total = 0;
        Double_t weight_sum = 0;
        for (Int_t i = 0; i < points.size(); ++i)
        {
            const Double_t d = DistanceToLine(points.at(i), linePoint, lineDir);
            d2total += hit_energy.at(i) * Power(d, 2);
            weight_sum += hit_energy.at(i);
        }
        radius = Sqrt(d2total / weight_sum);

        delete minimiser;
        points.clear();
    }
}

Double_t Fit3D::DistanceToLine(const TVector3& point, const TVector3& linePoint, const TVector3& lineDir)
{
    TVector3 vec = point - linePoint;
    return (vec - vec.Dot(lineDir.Unit()) * lineDir.Unit()).Mag();
}

Double_t Fit3D::SumDist2(const Double_t* par)
{
    TVector3 linePoint(par[0], par[1], par[2]);
    TVector3 lineDir(par[3], par[4], par[5]);
    lineDir = lineDir.Unit();

    Double_t sum = 0;
    for (Int_t i = 0; i < points.size(); ++i)
    {
        const Double_t d = DistanceToLine(points.at(i), linePoint, lineDir);
        sum += weights.at(i) * Power(d, 2);
    }
    return sum;
}
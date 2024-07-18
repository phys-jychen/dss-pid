#ifndef FIT3D_H
#define FIT3D_H

#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include "TROOT.h"
#include "TMath.h"
#include "TGraph2D.h"
#include "TF3.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TSystem.h"
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

using std::vector;
using namespace ROOT;
using namespace TMath;

class Fit3D
{
public:
    Fit3D(const vector<Double_t>& hit_x, const vector<Double_t>& hit_y, const vector<Double_t>& hit_z, const Int_t& nhits);

    ~Fit3D() = default;

    [[nodiscard]] Double_t GetRMSRadius() const { return radius; }

    [[maybe_unused]] [[nodiscard]] Double_t GetDirectionX() const { return direction_x; }

    [[maybe_unused]] [[nodiscard]] Double_t GetDirectionY() const { return direction_y; }

    [[maybe_unused]] [[nodiscard]] Double_t GetDirectionZ() const { return direction_z; }

private:
    Double_t radius = 0;

    Double_t direction_x = 0;

    Double_t direction_y = 0;

    Double_t direction_z = 0;

    static Double_t DistanceToLine(const TVector3& point, const TVector3& linePoint, const TVector3& lineDir);

    static Double_t SumDist2(const Double_t* par);
};

#endif

#pragma once

#include "MoveBase.h"
#include "System.h"
#include "StaticVals.h"

class MultiParticle : public MoveBase
{
public:
  MultiParticle(System &sys, StaticVals const& statV);

  virtual uint Prep(const double subDraw, const double movPerc);
  virtual void CalcEn();
  virtual uint Transform();
  virtual void Accept(const uint rejectState, const uint step);
  virtual double GetCoeff();
private:
  uint bPick;
  SystemPotential sysPotNew;
  XYZArray& atomForceRef;
  XYZArray atomForceNew;
  XYZArray& atomTorqueRef;
  XYZArray atomTorqueNew;
  XYZArray& molForceRef;
  XYZArray molForceNew;
  XYZArray& molTorqueRef;
  XYZArray molTorqueNew;
  XYZArray t_k;
  XYZArray r_k;
  Coordinates newMolsPos;
  COM newCOMs;
  vector<uint> moveType;
  const MoleculeLookup& molLookup;
  double w_new, w_old;
  double t_max, r_max;
  double lambda;
};

inline MultiParticle::MultiParticle(System &sys, StaticVals const &statV) :
  MoveBase(sys, statV),
  atomForceRef(sys.atomForceRef),
  atomTorqueRef(sys.atomTorqueRef),
  molForceRef(sys.molForceRef),
  molTorqueRef(sys.molTorqueRef),
  molLookup(sys.molLookup)
{
  atomForceNew.Init(sys.atomForceRef.Count());
  atomTorqueNew.Init(sys.atomTorqueRef.Count());
  molForceNew.Init(sys.molForceRef.Count());
  molTorqueNew.Init(sys.molTorqueRef.Count());
  t_k.Init(sys.com.Count());
  r_k.Init(sys.com.Count());
  newMolPos.Init(sys.coordinates.Count());
  newCOMs.Init(sys.com.Count());
  moveType.resize(sys.com.Count());
  lambda = 0.5;
}

inline uint MultiParticle::Prep(const double subDraw, const double movPerc)
{
  uint state = mv::fail_state::NO_FAIL;
#if ENSEMBLE == GCMC
  bPick = mv::BOX0;
#else
  prng.PickBox(bPick, subDraw, movPerc);
#endif
  // subPick = mv::GetMoveSubIndex(mv::MULTIPARTICLE, bPick);
  t_max = 0.05;
  r_max = 0.09 * 2 * M_PI;

  MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(bPick);
  MoleculeLookup::box_iterator end = molLook.BoxEnd(bPick);
  while(thisMol != end) {
    uint length = molRef.GetKind(*thisMol).NumAtoms();
    if(length==1)
      moveType[*thisMol] = 0;
    else
      moveType[*thisMol] = prng.RandInt(1);
  }
  coordCurrRef.CopyRange(newMolsPos, 0, 0, coordCurrRef.Count());
  comCurrRef.CopyRange(newCOMs, 0, 0, comCurrRef.Count());
  return state;
}

inline uint MultiParticle::Transform()
{
  // Based on the reference force decided whether to displace or rotate each
  // individual particle.

  // move particles according to force and torque and store them in the new pos
  return 0;
}

inline void MultiParticle::CalcEn() 
{
  // Calculate the new force and energy and we will compare that to the
  // reference values in Accept() function
  cellList.GridAll(boxDimRef, newMolsPos, molLookRef);

  sysPotNew = sysPotRef;
  sysPotNew = calcEnRef.BoxInter(sysPotNew, newMolPos, newCOMs, atomForceNew,
                                 molForceNew, atomTorqueNew, molTorqueNew,
                                 boxDimRef, bPick);
  return;
}

inline double MultiParticle::GetCoeff() const
{
  // calculate (w_new->old/w_old->new) and return it.
  uint length, start;
  XYZ lbf; // lambda * BETA * force
  XYZ lbt; // lambda * BETA * torque
  w_old = 1.0;
  w_new = 1.0;
  uint molNumber;
  MoleculeLookup::box_iterator thisMol = molLookup.BoxBegin(bPick);
  MoleculeLookup::box_iterator end = molLook.BoxEnd(bPick);

  while(thisMol != end) {
    molNumber = *thisMol;
    if(moveType[molNumber]) { // == 1 -> rotate
      lbt = molTorqueRef.Get(molNumber) * lambda * BETA;
      w_old *= lbt.x * exp(lbf.x * r_k.Get(molNumber]).x)/
        (2.0*sinh(lbt.x * r_max));
      w_old *= lbt.y * exp(lbf.y * r_k.Get(molNumber]).y)/
        (2.0*sinh(lbt.y * r_max));
      w_old *= lbt.z * exp(lbf.z * r_k.Get(molNumber]).z)/
        (2.0*sinh(lbt.z * r_max));

      lbt = molTorqueNew.Get(molNumber) * lambda * BETA;
      w_new *= lbt.x * exp(lbf.x * -1 * r_k.Get(molNumber]).x)/
        (2.0*sinh(lbt.x * r_max));
      w_new *= lbt.y * exp(lbf.y * -1 * r_k.Get(molNumber]).y)/
        (2.0*sinh(lbt.y * r_max));
      w_new *= lbt.z * exp(lbf.z * -1 * r_k.Get(molNumber]).z)/
        (2.0*sinh(lbt.z * r_max));
    }
    else { // displace
      lbf = molForceRef.Get(molNumber) * lambda * BETA;
      w_old *= lbf.x * exp(lbf.x * t_k.Get(molNumber]).x)/
        (2.0*sinh(lbf.x * t_max));
      w_old *= lbf.y * exp(lbf.y * t_k.Get(molNumber]).y)/
        (2.0*sinh(lbf.y * t_max));
      w_old *= lbf.z * exp(lbf.z * t_k.Get(molNumber]).z)/
        (2.0*sinh(lbf.z * t_max));

      lbf = molForceNew.Get(molNumber) * lambda * BETA;
      w_new *= lbf.x * exp(lbf.x * -1 * t_k.Get(molNumber]).x)/
        (2.0*sinh(lbf.x * t_max));
      w_new *= lbf.y * exp(lbf.y * -1 * t_k.Get(molNumber]).y)/
        (2.0*sinh(lbf.y * t_max));
      w_new *= lbf.z * exp(lbf.z * -1 * t_k.Get(molNumber]).z)/
        (2.0*sinh(lbf.z * t_max));

    }
  }
  return w_old/w_new;
}

inline void MultiParticle::Accept(const uint rejectState, const uint step)
{
  // Here we compare the values of reference and trial and decide whether to 
  // accept or reject the move
  double MPCoeff = GetCoeff();
  double uBoltz = exp(-BETA * (sysPotNew.Total() - sysPotRef.Total()));
  double accept = MPCoeff * uBoltz;
  bool result = (rejectState == mv::fail_state::NO_FAIL) && prng() < accept;
  if(result) {
    sysPotRef = sysPotNew;
    swap(coordCurrRef, newMolsPos);
    swap(comCurrRef, newCOMs);
    swap(molForceRef, molForceNew);
    swap(atomForceRef, atomForceNew);
    swap(molTorqueRef, molTorqueNew);
    swap(atomTorqueRef, atomTorqueNew);
  }
  else {
    cellList.GridAll(boxDimRef, coordCurrRef, molLookRef);
  }
  return;
}
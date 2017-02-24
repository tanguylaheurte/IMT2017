/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2003 Ferdinando Ametrano
 Copyright (C) 2001, 2002, 2003 Sadruddin Rejeb
 Copyright (C) 2005 StatPro Italia srl
 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/
 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.
 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file binomialtree.hpp
    \brief Binomial tree class
*/

#ifndef binomial_tree_hpp
#define binomial_tree_hpp

#include <ql/methods/lattices/tree.hpp>
#include <ql/instruments/dividendschedule.hpp>
#include <ql/stochasticprocess.hpp>
#include <iostream>

namespace QuantLib {

    //! Binomial tree base class
    /*! \ingroup lattices */
    template <class T>
    class BinomialTree_2 : public Tree<T> {
      public:
        enum Branches { branches = 2 };
        BinomialTree_2(const boost::shared_ptr<StochasticProcess1D>& process,
                       Time end,
                       Size steps)
        : Tree<T>(steps+1) {
            x0_ = process->x0(); std::cout 
<< "hello" << std::endl;
            dt_ = end/steps;
            driftPerStep_ = process->drift(0.0, x0_) * dt_;
        }
        Size size(Size i) const {;
            return i+3; // i+1->i+3 to get 3nods at t=0
        }
        Size descendant(Size, Size index, Size branch) const {
            return index + branch;
        }
      protected:
        Real x0_, driftPerStep_;
        Time dt_;
    };


    //! Base class for equal probabilities binomial tree
    /*! \ingroup lattices */
    template <class T>
    class EqualProbabilitiesBinomialTree_2 : public BinomialTree_2<T> {
      public:
        EqualProbabilitiesBinomialTree_2(
                        const boost::shared_ptr<StochasticProcess1D>& process,
                        Time end,
                        Size steps)
        : BinomialTree_2<T>(process, end, steps) {}
        Real underlying(Size i, Size index) const {
            BigInteger j = 2*BigInteger(index) - BigInteger(i); //We have to modify this to get the good amount of nods at each step
            // exploiting the forward value tree centering
            return this->x0_*std::exp(i*this->driftPerStep_ + j*this->up_);
        }
        Real probability(Size, Size, Size) const { return 0.5; }
      protected:
        Real up_;
    };


    //! Base class for equal jumps binomial tree
    /*! \ingroup lattices */
    template <class T>
    class EqualJumpsBinomialTree_2 : public BinomialTree_2<T> {
      public:
        EqualJumpsBinomialTree_2(
                        const boost::shared_ptr<StochasticProcess1D>& process,
                        Time end,
                        Size steps)
        : BinomialTree_2<T>(process, end, steps) {}
        Real underlying(Size i, Size index) const {
            BigInteger j = 2*BigInteger(index) - BigInteger(i); //We have to modify this to get the good amount of nods at each step
            // exploiting equal jump and the x0_ tree centering
            return this->x0_*std::exp(j*this->dx_);
        }
        Real probability(Size, Size, Size branch) const {
            return (branch == 1 ? pu_ : pd_);
        }
      protected:
        Real dx_, pu_, pd_;
    };


    //! Jarrow-Rudd (multiplicative) equal probabilities binomial tree
    /*! \ingroup lattices */
    class JarrowRudd_2 : public EqualProbabilitiesBinomialTree_2<JarrowRudd_2> {
      public:
        JarrowRudd_2(const boost::shared_ptr<StochasticProcess1D>&,
                     Time end,
                     Size steps,
                     Real strike);
    };


    //! Cox-Ross-Rubinstein (multiplicative) equal jumps binomial tree
    /*! \ingroup lattices */
    class CoxRossRubinstein_2
        : public EqualJumpsBinomialTree_2<CoxRossRubinstein_2> {
      public:
        CoxRossRubinstein_2(const boost::shared_ptr<StochasticProcess1D>&,
                            Time end,
                            Size steps,
                            Real strike);
    };


    //! Additive equal probabilities binomial tree
    /*! \ingroup lattices */
    class AdditiveEQPBinomialTree_2
        : public EqualProbabilitiesBinomialTree_2<AdditiveEQPBinomialTree_2> {
      public:
        AdditiveEQPBinomialTree_2(
                        const boost::shared_ptr<StochasticProcess1D>&,
                        Time end,
                        Size steps,
                        Real strike);
    };


    //! %Trigeorgis (additive equal jumps) binomial tree
    /*! \ingroup lattices */
    class Trigeorgis_2 : public EqualJumpsBinomialTree_2<Trigeorgis_2> {
      public:
        Trigeorgis_2(const boost::shared_ptr<StochasticProcess1D>&,
                     Time end,
                     Size steps,
                     Real strike);
    };


    //! %Tian tree: third moment matching, multiplicative approach
    /*! \ingroup lattices */
    class Tian_2 : public BinomialTree_2<Tian_2> {
      public:
        Tian_2(const boost::shared_ptr<StochasticProcess1D>&,
               Time end,
               Size steps,
               Real strike);
        Real underlying(Size i, Size index) const {
            return x0_ * std::pow(down_, Real(BigInteger(i)-BigInteger(index)))
                       * std::pow(up_, Real(index));
        };
        Real probability(Size, Size, Size branch) const {
            return (branch == 1 ? pu_ : pd_);
        }
      protected:
        Real up_, down_, pu_, pd_;
    };

    //! Leisen & Reimer tree: multiplicative approach
    /*! \ingroup lattices */
    class LeisenReimer_2 : public BinomialTree_2<LeisenReimer_2> {
      public:
        LeisenReimer_2(const boost::shared_ptr<StochasticProcess1D>&,
                       Time end,
                       Size steps,
                       Real strike);
        Real underlying(Size i, Size index) const {
            return x0_ * std::pow(down_, Real(BigInteger(i)-BigInteger(index)))
                       * std::pow(up_, Real(index));
        }
        Real probability(Size, Size, Size branch) const {
            return (branch == 1 ? pu_ : pd_);
        }
      protected:
        Real up_, down_, pu_, pd_;
    };


     class Joshi4_2 : public BinomialTree_2<Joshi4_2> {
      public:
        Joshi4_2(const boost::shared_ptr<StochasticProcess1D>&,
                 Time end,
                 Size steps,
                 Real strike);
        Real underlying(Size i, Size index) const {
            return x0_ * std::pow(down_, Real(BigInteger(i)-BigInteger(index)))
                       * std::pow(up_, Real(index));
        }
        Real probability(Size, Size, Size branch) const {
            return (branch == 1 ? pu_ : pd_);
        }
      protected:
        Real computeUpProb(Real k, Real dj) const;
        Real up_, down_, pu_, pd_;
    };

}


#endif

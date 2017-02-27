#include "binomialtree.hpp"
#include "binomialengine.hpp"
#include <ql/stochasticprocess.hpp>

#include <ql/pricingengines/vanilla/all.hpp>
#include <iostream>
#include <ql/time/calendars/target.hpp>
#include <boost/timer.hpp>
#include <ctime> 

using namespace QuantLib;

int main() {

    try {
		 
        std::cout << std::endl;
		Size timeSteps = 801;
        // set up dates
        Calendar calendar = TARGET();
        Date todaysDate(15, May, 1998);
        Date settlementDate(18, May, 1998);
        Settings::instance().evaluationDate() = todaysDate;

        // our options
        Option::Type type(Option::Put);
        Real underlying = 36;
        Real strike = 40;
        Spread dividendYield = 0.00;
        Rate riskFreeRate = 0.06;
        Volatility volatility = 0.20;
        Date maturity(17, May, 2002);
        DayCounter dayCounter = Actual365Fixed();

        std::cout << "Option type = "  << type << std::endl;
        std::cout << "Maturity = "        << maturity << std::endl;
        std::cout << "Underlying price = "        << underlying << std::endl;
        std::cout << "Strike = "                  << strike << std::endl;

        std::cout << std::endl;
        std::string method;
        std::cout << std::endl ;

        // write column headings
        Size widths[] = { 35, 14, 14, 14 };
        std::cout << std::setw(widths[0]) << std::left << "Method"
                  << std::setw(widths[1]) << std::left << "European"
                  << std::endl;

        std::vector<Date> exerciseDates;
        for (Integer i=1; i<=4; i++)
            exerciseDates.push_back(settlementDate + 3*i*Months);

        boost::shared_ptr<Exercise> europeanExercise(
                                         new EuropeanExercise(maturity));
		
		Handle<Quote> underlyingH(
            boost::shared_ptr<Quote>(new SimpleQuote(underlying)));
		// bootstrap the yield/dividend/vol curves
        Handle<YieldTermStructure> flatTermStructure(
            boost::shared_ptr<YieldTermStructure>(
                new FlatForward(settlementDate, riskFreeRate, dayCounter)));
        Handle<YieldTermStructure> flatDividendTS(
            boost::shared_ptr<YieldTermStructure>(
                new FlatForward(settlementDate, dividendYield, dayCounter)));
        Handle<BlackVolTermStructure> flatVolTS(
            boost::shared_ptr<BlackVolTermStructure>(
                new BlackConstantVol(settlementDate, calendar, volatility,
                                     dayCounter)));
        boost::shared_ptr<StrikedTypePayoff> payoff(
                                        new PlainVanillaPayoff(type, strike));
        boost::shared_ptr<BlackScholesMertonProcess> bsmProcess(new BlackScholesMertonProcess(underlyingH, flatDividendTS,
                                               flatTermStructure, flatVolTS));

        // options
        VanillaOption europeanOption(payoff, europeanExercise);
		//B&S=3.844308
		//OldBinomial=3.843504
		//newBinomial=3.843431
 
        // Analytic formulas:

        // Black-Scholes for European
        method = "Black-Scholes";
        europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(new AnalyticEuropeanEngine(bsmProcess)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::endl;
		
		Real deltaBS=europeanOption.delta();
 		Real gammaBS=europeanOption.gamma();

		std::cout << "delta via BS: " << deltaBS << std::endl;
		std::cout << "Gamma via BS: " << gammaBS << std::endl << std::endl;
		
		//start chrono		
		const clock_t begin_time = clock();
		


		method = "Binomial Cox-Ross-Rubinstein_2";
        europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
                      new BinomialVanillaEngine_2<CoxRossRubinstein_2>(bsmProcess,
                                                                   timeSteps)));
		
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::endl;
		
		Real deltaBinomial=europeanOption.delta();
 		Real gammaBinomial=europeanOption.gamma();
		Real diffDelta=deltaBinomial-deltaBS;
 		Real diffgamma=gammaBinomial-gammaBS;

		std::cout << "delta via BinomialCRR: " << deltaBinomial << std::endl;	

		std::cout << "Gamma via BinomialCRR: " << gammaBinomial << std::endl << std::endl;
		
		std::cout << "deltaBinomialCRR - deltaBS = " << diffDelta << std::endl;
		std::cout << "gammaBinomialCRR - gammaBS = " << diffgamma << std::endl << std::endl;
	
		std::cout << "Time spent to build the tree in sec : " <<  float( clock () - begin_time) /  CLOCKS_PER_SEC << std::endl << std::endl;

		//start chrono		
		const clock_t begin_time2 = clock();

		method = "Joshi_2";
        europeanOption.setPricingEngine(boost::shared_ptr<PricingEngine>(
                      new BinomialVanillaEngine_2<Joshi4_2>(bsmProcess,
                                                                   timeSteps)));
		
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::endl;
		
		Real deltaBinomialJ=europeanOption.delta();
 		Real gammaBinomialJ=europeanOption.gamma();
		diffDelta=deltaBinomialJ-deltaBS;
 		diffgamma=gammaBinomialJ-gammaBS;

		std::cout << "delta via BinomialJ: " << deltaBinomialJ << std::endl;	

		std::cout << "Gamma via BinomialJ: " << gammaBinomialJ << std::endl << std::endl;
		
		std::cout << "deltaBinomialJ - deltaBS = " << diffDelta << std::endl ;
		std::cout << "gammaBinomialJ - gammaBS = " << diffgamma << std::endl << std::endl;
	
		std::cout << "Time spent to build the tree in sec : " <<  float( clock () - begin_time2) /  CLOCKS_PER_SEC << std::endl << std::endl;
        return 0;

    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "unknown error" << std::endl;
        return 1;
    }
}

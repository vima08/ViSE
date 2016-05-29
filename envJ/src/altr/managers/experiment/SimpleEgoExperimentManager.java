/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package altr.managers.experiment;

import altr.Environment;
import static altr.Runner.writer;
import altr.distributions.Distribution;
import altr.entity.Group;
import altr.entity.Offer;
import altr.entity.Person;
import altr.experiment.Experiment;
import altr.managers.Analyzer;
import altr.managers.PersonManager;
import altr.strategies.SimpleEgoisticStrategy;
import altr.strategies.SimpleGroupStrategy;
import java.util.ArrayList;
import java.util.Collection;

/**
 *
 * @author Vitaly
 */
public class SimpleEgoExperimentManager extends ExperimentManager {    
    
    Collection<Person> others;    
    double[] othersAvgMoney;    
    double[] acceptanceCounters;
    double alpha; // acceptance
    SimpleGroupStrategy groupStrategy;

    public SimpleEgoExperimentManager(Experiment experiment, Environment env, Person egoist, int egoSize, Double alpha) throws CloneNotSupportedException {
        super(experiment, env);                
        this.alpha = alpha;
        othersAvgMoney = new double[stepNumber];         
        acceptanceCounters = new double[stepNumber]; 
        
        egoist.setStrategy(new SimpleEgoisticStrategy());
        
        this.others = PersonManager.clonePerson(egoist, egoSize);
        
        pM.getPeople().addAll(others);
    }
    
    @Override
    protected void accept(Boolean isAccepted, Collection<Offer> offers) {
        if (!isAccepted) return;
        ArrayList<Person> p = (ArrayList)pM.getPeople();
        ArrayList<Offer> o = (ArrayList)offers;
        for(int i =0; i < pM.getPeople().size(); i++) {
            Double money = p.get(i).getMoney();
            p.get(i).setMoney(money + o.get(i).getAmmount());
        }
    }

    @Override
    protected void results() {      
//        writer.comment("test");
//        writer.value("Average capital").newLine();
//        for(int i = 0; i < stepNumber; i++){
//            writer.value(Double.toString(groupAvgMoney[i]/iterationNumber)).value(Double.toString(othersAvgMoney[i]/iterationNumber)).
//                    value(Double.toString(avgMoney[i]/iterationNumber)).value(Double.toString(acceptanceCounters[i]/iterationNumber)).newLine();
//        }
        int last = stepNumber-1;
        writer.value(Double.toString(othersAvgMoney[last]/iterationNumber)).
            value(Double.toString(acceptanceCounters[last]/iterationNumber)).
            newLine();
        System.out.println("Finished!");
    }

    @Override
    protected void analysis() {        
        double avg1 = Analyzer.getAverageMoney(others);
        othersAvgMoney[step] += avg1;
        acceptanceCounters[step] += acceptanceCounter;
    }

    @Override
    protected Collection<Offer> generateOffers(Distribution dist, Collection<Person> people) {
        Collection<Offer> offers = new ArrayList();
        for(Person p: people) {
            Offer offer = new Offer(dist.getValue(), p.getId());
            offers.add(offer);
            //System.out.println(offer);
        }
        return offers;
    }

    @Override
    protected boolean isAccepted(Collection<Offer> offers) {
        double votes = 0;
        long number = pM.getPeople().size();
        for(Person p: pM.getPeople()) {
            if (p.getStrategy().vote(offers, pM.getPeople(),  p.getId())) votes++;
        }
        double percentage = votes / number;
        return (percentage >= alpha);
    }

    @Override
    protected void postprocessing() {        
        for(Person p: pM.getPeople()){
            p.setActive(p.getMoney() > 0);
        }
    }
    
}

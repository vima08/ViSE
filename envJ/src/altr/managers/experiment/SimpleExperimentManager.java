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
public class SimpleExperimentManager extends ExperimentManager {    
       
    Collection<Person> group;
    Collection<Person> others;
    double[] groupAvgMoney;
    double[] othersAvgMoney;
    double[] avgMoney;
    double[] acceptanceCounters;
    int groupSize;
    SimpleGroupStrategy groupStrategy;

    public SimpleExperimentManager(Experiment experiment, Environment env, Person groupMan, int groupSize, Person egoist, int egoSize) throws CloneNotSupportedException {
        super(experiment, env);
        this.groupSize = groupSize;
        groupAvgMoney = new double[stepNumber]; 
        othersAvgMoney = new double[stepNumber]; 
        avgMoney = new double[stepNumber]; 
        acceptanceCounters = new double[stepNumber]; 
        groupStrategy = new SimpleGroupStrategy(gM, 0.0);
        groupMan.setStrategy(groupStrategy);
        egoist.setStrategy(new SimpleEgoisticStrategy());
        this.group = PersonManager.clonePerson(groupMan, groupSize);
        this.others = PersonManager.clonePerson(egoist, egoSize);

        Group g = new Group("test group", "group1");
        gM.addGroup(g);
        gM.addPeopleToGroup(group, g.getId());
        groupStrategy.init();
        pM.getPeople().addAll(group);
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
        writer.value(Double.toString(groupAvgMoney[last]/iterationNumber)).value(Double.toString(othersAvgMoney[last]/iterationNumber)).
            value(Double.toString(avgMoney[last]/iterationNumber)).value(Double.toString(acceptanceCounters[last]/iterationNumber)).
            newLine();
        System.out.println("Finished!");
    }

    @Override
    protected void analysis() {
        double avg = Analyzer.getAverageMoney(group);
        double avg1 = Analyzer.getAverageMoney(others);
        double avg2 = Analyzer.getAverageMoney(pM.getPeople());
        groupAvgMoney[step] += avg;
        othersAvgMoney[step] += avg1;
        avgMoney[step] += avg2;  
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
        return (percentage >= 0.5);
    }

    @Override
    protected void postprocessing() {
        groupStrategy.reset();
        for(Person p: pM.getPeople()){
            p.setActive(p.getMoney() > 0);
        }
    }
    
}

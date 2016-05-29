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
import altr.strategies.hierarchical.HierarchicalBossStrategy;
import altr.strategies.hierarchical.HierarchicalEmployeeStrategy;
import altr.strategies.hierarchical.HierarchicalManagerStrategy;
import java.util.ArrayList;
import java.util.Collection;

/**
 *
 * @author Vitaly
 */
public class HierarchicalExperimentM extends ExperimentManager {    
    
    long m = 10;
    long n = 3;
       
    Collection<Person> boss;    
    Collection<Person> managers;
    Collection<Person> employees;
    double[] bossAvgMoney;
    double[] managersAvgMoney;
    double[] employeeAvgMoney;
    double[] avgMoney;
    double[] acceptanceCounters;    

    public HierarchicalExperimentM(Experiment experiment, Environment env, double alfa) throws CloneNotSupportedException {
        super(experiment, env);
        Person boss = new Person(10, true, 1, 0, 0, "boss", new HierarchicalBossStrategy(alfa));
        Person manager = new Person(10, true, 1, 0, 0, "manager", new HierarchicalManagerStrategy(alfa));
        Person employee = new Person(10, true, 1, 0, 0, "employee", new HierarchicalEmployeeStrategy(alfa));
        this.boss = PersonManager.clonePerson(boss, 1);
        this.managers = PersonManager.clonePerson(manager, 3);
        this.employees = PersonManager.clonePerson(employee, 9);
        
        pM.getPeople().addAll(this.boss);
        pM.getPeople().addAll(managers);
        pM.getPeople().addAll(employees);
        
        ArrayList<Person> people =  (ArrayList<Person>)pM.getPeople();
        for(int i = 0; i < 13; i++){
            people.get(i).setId(i);
        }
                
        bossAvgMoney = new double[stepNumber]; 
        employeeAvgMoney = new double[stepNumber]; 
        managersAvgMoney = new double[stepNumber];
        avgMoney = new double[stepNumber]; 
        acceptanceCounters = new double[stepNumber];                 
        
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
        writer.value(Double.toString(bossAvgMoney[last]/iterationNumber)).value(Double.toString(managersAvgMoney[last]/iterationNumber))
                .value(Double.toString(employeeAvgMoney[last]/iterationNumber)).
            value(Double.toString(avgMoney[last]/iterationNumber)).value(Double.toString(acceptanceCounters[last]/iterationNumber)).
            newLine();
        System.out.println("Finished!");
    }

    @Override
    protected void analysis() {
        double avg = Analyzer.getAverageMoney(boss);
        double avg1 = Analyzer.getAverageMoney(managers);        
        double avg2 = Analyzer.getAverageMoney(employees);
        double avg3 = Analyzer.getAverageMoney(pM.getPeople());
        bossAvgMoney[step] += avg;
        managersAvgMoney[step] += avg1;
        employeeAvgMoney[step] += avg2;
        avgMoney[step] += avg3;  
        acceptanceCounters[step] += acceptanceCounter;
    }

    @Override
    protected Collection<Offer> generateOffers(Distribution dist, Collection<Person> people) {
        Collection<Offer> offers = new ArrayList();
        for(Person p: people) {
            Offer offer = new Offer(dist.getValue(), p.getId());
            if ("boss".equals(p.getName())) {
                offer.setAmmount(offer.getAmmount() * 9);
            } else if ("manager".equals(p.getName())) {
                offer.setAmmount(offer.getAmmount() * 3);
            }
            offers.add(offer);
            //System.out.println(offer);
        }
        return offers;
    }

    @Override
    protected boolean isAccepted(Collection<Offer> offers) {
        double votes = 0;
        long number = this.employees.size() * 3;
        for(Person p: pM.getPeople()) {
            if (p.getStrategy().vote(offers, pM.getPeople(),  p.getId())) {
                if ("boss".equals(p.getName())) {
                    votes += 9;
                } else if ("manager".equals(p.getName())) {
                    votes += 3;
                } else {
                    votes += 1;
                }
            }
        }
        double percentage = votes / number;
        return (percentage >= 0.5);
    }

    @Override
    protected void postprocessing() {        
        for(Person p: pM.getPeople()){
            p.setActive(p.getMoney() > 0);
        }
    }
    
}

#include "gravitacek2/setup.hpp"
#include "gravitacek2/odesolver/ode.hpp"
#include "gravitacek2/odesolver/stepper.hpp"
#include "gravitacek2/odesolver/stepcontroller.hpp"
#include "gravitacek2/odesolver/event.hpp"

#include <vector>
#include <string>

namespace gr2
{
    class Integrator
    {
    protected:
        // ========== Components of itnegrator ========== 
        ODE *ode;                       // ordinary differential equation
        Stepper *stepper;               // stepper used for integrating
        StepController *stepcontroller; // step controller

        // ========== Events ========== 
        std::vector<Event*> events_data;
        std::vector<Event*> events_modifying;
        std::vector<Event*> events_terminal;

        // ========== Values of events ========== 
        real* events_data_values;
        real* events_modifying_values;
        real* events_terminal_values;

        // ========== Number of events ========== 
        int number_of_events_data;
        int number_of_events_modifying;
        int number_of_events_terminal;

        // ========== Scalar values ========== 
        real t;                         // current value of time
        real h;                         // current value of time step
        real h2;
        real h3;

        // ========== Arrays of values ========== 
        real *yt;                       // array for storing current value of y
        real *yt2;                      // array for trying next step
        real *yt3;                      // array for trying events
        real *dydt;                     // array for calculating \f$\frac{\mathrm{d} \vec{y}}{\mathrm{d}}\f$
        real *dydt2;                    // 
        real *dydt3;
        real *err;
        real *err2;
        real *err3;

        void basic_setup();

        /**
         * @brief Ininialize stepper based on its name.
         * 
         * @param stepper_name Name of stepper
         */
        void init_stepper(const std::string& stepper_name);
        bool solve_event(Event *event, real& previous_value_of_event);

    public:
        /**
         * @brief Construct a new Integrator object.
         * 
         * This constructor assumes constant steps in time.
         * 
         * @param ode set of ordinary differential equations to be solved
         * @param stepper_name name of stepper
         */
        Integrator(ODE &ode, const std::string& stepper_name);

        /**
         * @brief Construct a new Integrator object.
         * 
         * This constructor assumes StandardStepController.
         * 
         * @param ode set of ordinary differential equations to be solved
         * @param stepper_name name of stepper
         * @param eps_abs expected absolute error
         * @param eps_rel expected relative error
         * @param a_y coefficient 
         * @param a_dydt coefficient
         */
        Integrator(ODE &ode, const std::string& stepper_name, const real &atol, const real &rtol);

        /**
         * @brief Destroy the Integrator object.
         * 
         */
        ~Integrator();

        /**
         * @brief Add event to simulation.
         * 
         * @param event 
         */
        void add_event(Event* event);

        /**
         * @brief Integrate ordinary differential equation
         * 
         * @param y_start 
         * @param t_start 
         * @param t_end 
         * @param h_start 
         */
        void integrate(const real y_start[], const real &t_start, const real &t_end, const real &h_start);
    };
}
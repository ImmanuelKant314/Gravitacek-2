#include "gravitacek2/integrator/event.hpp"

namespace gr2
{
    Event::Event(const EventType &type, const bool &terminal):type(type), terminal(terminal)
    {
    }

    EventType Event::get_type() const 
    {
        return this->type;
    }

    bool Event::get_terminal() const
    {
        return this->terminal;
    }
}
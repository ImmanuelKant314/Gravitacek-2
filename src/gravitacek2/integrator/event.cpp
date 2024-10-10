#include "gravitacek2/integrator/event.hpp"

namespace gr2
{
    Event::Event(const EventType& type)
    {
        this->type = type;
    }

    EventType Event::get_type() const 
    {
        return this->type;
    }
}
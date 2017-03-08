# Double pendulum with a single collision pt
Hybrid system: continuous integration with discrete collision events

Run MAIN.m

Integration switches between 3 equations of motion depending on the
contact configuration:
1. Free motion
2. Collision with link 1 -- pt contact
3. Collision with link 2 -- sliding contact with friction

Contact is broken when normal force goes negative OR link 2 slides off
the contact pt.
Collisions can be either plastic or elastic.
Only 1 contact pt at a time

Deriver.m symbolically derives the equations of motion, the contact
conditions, collision resolutions, etc, and automatically writes them to
files. MAIN.m, Deriver, RHS.m, and eventFun.m are the only hand-written
functions.

Matthew Sheen, 2017

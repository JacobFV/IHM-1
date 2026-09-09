"""Uncalibrated rigid contact materials and passive sphere/plane impact solver."""
import numpy as np


BOUNCY_BALL = {'name': 'engineering-rubber', 'restitution': .75, 'friction': .35,
               'bounce_threshold_m_s': .05}


def sphere_plane_step(x, v, omega, dt, acceleration, radius, mass, axis, plane,
                      restitution=.75, friction=.35, bounce_threshold=.05):
    """Advance constant acceleration with time-of-impact splitting, in place.

    Normal restitution and clamped Coulomb tangential impulses never increase
    kinetic energy. Low-speed contact rests instead of numerically chattering.
    Position repair is reported separately and never converted into velocity.
    """
    inertia=.4*mass*radius**2; normal=np.eye(3)[axis];lever=-radius*normal
    a=np.array(acceleration,float);remaining=float(dt);loss=0.;impulse=np.zeros(3)
    correction=max(0.,plane+radius-x[axis]);x[axis]+=correction
    for _ in range(32):
        height=max(0.,x[axis]-plane-radius)
        if height<1e-12 and v[axis]>=0 and v[axis]<bounce_threshold and a[axis]<=0:
            # A low-energy contact is held by an explicit ideal plane reaction.
            loss+=.5*mass*v[axis]**2;impulse[axis]-=mass*v[axis]
            v[axis]=0.;impulse[axis]-=mass*a[axis]*remaining;a[axis]=0.
            # The support's Coulomb budget also damps sliding at rest.
            relative=v+np.cross(omega,lever);tangent=relative-relative[axis]*normal
            speed=np.linalg.norm(tangent);budget=friction*max(0.,impulse[axis])
            if speed:
                jt=-min(speed/(1/mass+radius**2/inertia),budget)*tangent/speed
                before=.5*mass*(v@v)+.5*inertia*(omega@omega)
                v+=jt/mass;omega+=np.cross(lever,jt)/inertia;impulse+=jt
                loss+=before-(.5*mass*(v@v)+.5*inertia*(omega@omega))
            x+=v*remaining+.5*a*remaining**2;v+=a*remaining
            return loss,impulse,correction
        if height<1e-12 and v[axis]<0:hit=0.
        elif abs(a[axis])<1e-15:hit=-height/v[axis] if v[axis]<0 else np.inf
        else:
            discriminant=v[axis]**2-2*a[axis]*height
            roots=[] if discriminant<0 else [(-v[axis]+s*np.sqrt(discriminant))/a[axis] for s in (-1,1)]
            hit=min((t for t in roots if t>1e-12 and v[axis]+a[axis]*t<0),default=np.inf)
        elapsed=min(hit,remaining)
        x+=v*elapsed+.5*a*elapsed**2;v+=a*elapsed;remaining-=elapsed
        if hit>elapsed:return loss,impulse,correction
        x[axis]=plane+radius
        relative=v+np.cross(omega,lever);vn=float(relative[axis]);e=restitution if -vn>=bounce_threshold else 0.
        jn=max(0.,-(1+e)*vn*mass);tangent=relative-vn*normal;speed=np.linalg.norm(tangent)
        jt=-min(speed/(1/mass+radius**2/inertia),friction*jn)*tangent/speed if speed else np.zeros(3)
        j=jn*normal+jt;before=.5*mass*(v@v)+.5*inertia*(omega@omega)
        v+=j/mass;omega+=np.cross(lever,j)/inertia;impulse+=j
        loss+=before-(.5*mass*(v@v)+.5*inertia*(omega@omega))
        if remaining<=1e-15:return loss,impulse,correction
    raise RuntimeError('Sphere contact exceeded impact iteration budget')

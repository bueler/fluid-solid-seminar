---
title: Math 692 Graduate Seminar on Fluids and Solids
---

This is the website for **Math 692 Graduate Seminar** on **Fluids and Solids: Mathematics and Numerics** in Spring 2025!  We are in the [Dept. of Mathematics and Statistics](http://www.uaf.edu/dms/) at the [University of Alaska Fairbanks](http://www.uaf.edu/).

The official course schedule says "Fluids: Mathematics and Numerics", but I changed the title so elastic solids are in scope.  (I am too lazy to make the change official.)

Course details:
  * Organizer/Instructor: [Ed Bueler](http://bueler.github.io/), [elbueler@alaska.edu](mailto:elbueler@alaska.edu).
  * Time and place: Thursdays 3:40-4:40pm, Chapman 206.
  * In-person is preferred if you are on campus!  I will set up Zoom only if it makes sense.
  * Credits (CRN 35130): 1.0, but **non-credit attendance is welcomed and encouraged**.

[A schedule of topics appears at the bottom.](#schedule)  This schedule is subject to change!

## guiding principles

1. show up
2. talk about fluids and solids and stuff
3. try running stuff that does stuff

## content

My idea for this seminar is that there are several math graduate students already doing fluids or elasticity research.  Broader knowledge, with some kind of **focus on the mathematical structures underlying continuum problems**, and on **finite element/volume formulation of the equations**, is helpful to all of us.

I figure I might start the seminar with some kind of continuum mechanics introduction.  However, general continuum mechanics quickly gets ... boring?  In any case, there are many on campus who are expert in various specific fluids/solids problems, and **they are all invited**.  Anyone who is just curious is also invited.  Ideally the experts can start their presentations from a somewhat-common starting point, before getting to their particular models.

Tsunamis, linear elasticity, seismology, vulcanology, glaciers, magneto-hydrodynamics, thermodynamics, and water balloons are all in scope.  So are talks about the underlying mathematics, for instance Stokes theorem, differential forms, and Sobolev spaces.

Also, it is now pretty easy to do numerical simulations of fluids and solids using the finite element method, for instance using [Firedrake](https://www.firedrakeproject.org/index.html).  That open source Python-based library was the subject of last year's [finite element seminar](https://bueler.github.io/fe-seminar/), so we have examples to work from.

## <a id="schedule"></a> topics (schedule)

| Date   | Speaker            | Topic                                          | Links        |
|--------|--------------------|------------------------------------------------|--------------|
| 16 Jan | Ed Bueler          | continuum mechanics, toward Navier-Stokes      | [slides](slides/bueler16jan.pdf) <br> [recording](https://alaska.zoom.us/rec/share/UBUDwv4neSnh6j6_DgyAZf_ym1o8pLba5AFBeLxbvZONa3VuoCeXQ0nguA-u3Js1.3BkVIVgSWyBGjvN9?pwd=cZJ2KPRKJmn8sxqJwiaCOh8gWFinx9m5) |
| 23 Jan | Ed Bueler          | Navier-Stokes for incompressible, viscous fluids | [slides](slides/bueler23jan.pdf) <br>[recording](https://alaska.zoom.us/rec/share/1Lx_GvcoeddAthqP_uRJnU1yLzFjYKUT5tukvHSc7PPIh8khrNhnqrZQ8-J7kI3z.OY75--IaccrJpUHb?pwd=6_kCPJzT4hKuzsrIAllNPVwJ2c4G1p5Q) <br> [navierstokes.py](py/bueler/navierstokes.py) <br> [cavity.py](py/bueler/cavity.py) |
| 30 Jan | Nick Harrison      | irrotational and potential flows | [recording](https://alaska.zoom.us/rec/share/DkJg5URpezBJu5wk-BvhERRsoySJQ1NIduQBHBdTiAHpRWecBwbSuuQV7Eww-kc5.e2i8PL_9KHi6Q1kR?pwd=FBlRZUo9HIFLt_3TH-IUH8S9LTlbmjQz) |
|  6 Feb | Martin Truffer     | constitutive relations or how to deal with stress | [recording](https://alaska.zoom.us/rec/share/5QYFXLJ6etSPOD3Zy0QZYmfMEeginWf75orMK7miHAhIeGhLDoEwiRlEhgio_ez5.wUE4X_XPkVBRpP9y?pwd=_w4qjBA2_ufI55o23vcFpAMTlmGsW-K4) <br> [video](https://www.youtube.com/watch?v=UEB39-jlmdw) |
| 13 Feb | Ed Bueler          | Firedrake for Navier-Stokes | [slides](slides/bueler13feb.pdf) <br> [recording](https://alaska.zoom.us/rec/share/4JzHa_7LjdWHuWzICj2AgrbCWekjKeZDGSse8LfkMPppCazb9X0PEBfbt1jqoT0.8HmOqUGIOT2_Y19s?pwd=Zql8r96yMhBmdp9d-u3Xts1ZNVM8nmND) <br> [navierstokes.py](py/bueler/navierstokes.py) <br> [cavity.py](py/bueler/cavity.py) <br> [cylinder.geo](py/bueler/cylinder.geo) <br> [cylinder.py](py/bueler/cylinder.py) |
| 20 Feb | Amy Jenson         | conservation of energy in glaciers | [recording](https://alaska.zoom.us/rec/share/ul4JlLS9_P6Co_uTJZiuTG5IZbKVaR74wKlHyn9F1oUpSOJPWZQTGc6M7gsmEoDp.ZB8zT46AQCgTBMk8?pwd=iN0uTpEXNPRbbHi0p7uBoYDDS7-ybzW2) |
| 27 Feb | Austin Smith       | introduction to magneto-hydro-dynamics | [recording](https://alaska.zoom.us/rec/share/d0vXlwjnTC5BKgfZVkme2qk8h-Jw5FB4RXnU59a4nwS6onpkMUIkJY5jLahq51gS.L02EJ6kAxt0EJA_x?pwd=MzugvNa3rvGRbU9Z4miNNkRzQfQ3PmMS) |
|  6 Mar |  |  |  |
| 13 Mar | _Spring Break!_    |  |  |
| 20 Mar |  |  |  |
| 27 Mar |  |  |  |
|  3 Apr |  |  |  |
| 10 Apr |  |  |  |
| 17 Apr | Stefano Fochesatto | adaptive meshing for fluids |  |
| 24 Apr |  |  |  |

<!--
30 Jan | Ed Bueler | reference configuration, linear elasticity |
-->

### Github repo with this website, codes, and PDF slides

The Github repo is [github.com/bueler/fluid-solid-seminar](https://github.com/bueler/fluid-solid-seminar).  It contains this Jekyll website, and a `slides/` folder, and codes in `py/` etc folders.

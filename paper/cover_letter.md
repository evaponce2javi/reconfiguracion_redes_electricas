# Carta al editor (borrador)

> Borrador para EPSR. Completar los corchetes y ajustar el tono. Verificar en el sistema de envío si se pide sugerir
> revisores, y dejar fuera a cualquiera con conflicto de interés.

---

Dear Editor,

We submit the manuscript "Coverage-guaranteed probabilistic bounds on voltages and branch currents from sparse
affine forms and the range of the second-order term" for consideration as a research article in *Electric Power
Systems Research*.

Linearized probabilistic power flow underlies much of today's uncertainty analysis and chance-constrained operation.
Its accuracy is usually assessed through moments or density curves. We check it where operational decisions are
taken: in the extreme tails, against nonlinear AC power-flow samples, on four test systems of up to 2869 buses.

The manuscript makes three points:

1. **Linearized bounds fail in the tails, and moments do not fix them.** Cumulant and point-estimate methods leave up
   to 12% of the AC samples outside bus-voltage bounds designed for 0.27%, and up to 100% outside branch-current
   bounds. Correcting the linearization with the mean and variance of the second-order term, which is what published
   second-order methods do, is still not enough. The reason is that the linearization error is strongly skewed and
   heavy-tailed: for 98.9% of the buses of PEGASE-2869 its long tail points toward lower voltages.
2. **Using the range of the exact second-order term instead of its moments restores coverage.** We combine it with a
   distribution-free concentration radius on the sparse affine form of each output. The range is sampled from a
   closed-form quadratic model, so no power flow is solved, and an order-statistics argument makes the coverage level
   provable for the second-order model.
3. **The construction is validated well beyond the design point.** The validation covers stressed loading up to 95%
   of the loadability limit, correlated sources and non-uniform input distributions. It is applied to bus voltages,
   branch voltage drops and branch currents, and the limits of the approach are measured and reported. On
   PEGASE-2869 the bounds for all buses, drops and branch ends take about 10 s.

Second-order probabilistic load flow is not new: we build on, and cite, the quadratic and second-order methods
published in this area. What the manuscript adds is the evidence that their moment-based use does not reach the
coverage required at small violation levels, a bound that does, and a provable level for it.

The manuscript is original, has not been published, and is not under consideration elsewhere. All authors have
approved it and declare no competing interests. The code and data that reproduce every table and figure will be made
available [repository / DOI].

Sincerely,

[Name], on behalf of the authors
Escuela de Ingeniería Informática, Pontificia Universidad Católica de Valparaíso, Chile
[e-mail]

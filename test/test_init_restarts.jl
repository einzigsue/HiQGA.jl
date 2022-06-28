## 1D functions
Random.seed!(200)
x = readdlm("testing_data/testfunc.txt", ',', Float64, '\n')[:,1]
y = readdlm("testing_data/testfunc.txt", ',', Float64, '\n')[:,2]
σ, fractrain = 0.275, 0.5
ntrain = round(Int, fractrain*length(y))
ynoisy = similar(y) .+ NaN
linidx = randperm(length(y))[1:ntrain]
ynoisy[linidx] = y[linidx] + σ*randn(ntrain)
line = transD_GP.Line(ynoisy;useML=false, σ=σ)

## make options for the multichannel lengthscale GP
log10bounds = [-1.2 -0.8]
δlog10λ = 0.05
nminlog10λ, nmaxlog10λ = 2, 30
pnorm = 2.
Klog10λ = transD_GP.GP.Mat32()
λlog10λ = [0.02abs(diff([extrema(x)...])[1])]
demean = false
sdev_poslog10λ = [0.05abs(diff([extrema(x)...])[1])]
sdev_proplog10λ = 0.05*diff(log10bounds, dims=2)[:]
xall = permutedims(collect(x))
xbounds = permutedims([extrema(x)...])
## Initialize a lengthscale model using these options
Random.seed!(12)
optlog10λ = transD_GP.OptionsStat(nmin = nminlog10λ,
                        nmax = nmaxlog10λ,
                        xbounds = xbounds,
                        fbounds = log10bounds,
                        xall = xall,
                        λ = λlog10λ,
                        δ = δlog10λ,
                        demean = demean,
                        sdev_prop = sdev_proplog10λ,
                        sdev_pos = sdev_poslog10λ,
                        needλ²fromlog = true,
                        updatenonstat = true,
                        pnorm = pnorm,
                        quasimultid = false,
                        K = Klog10λ,
                        save_freq = 20,
                        timesλ = 4,
                        peskycholesky = true
                        )
## make options for the nonstationary actual properties GP
nmin, nmax = 2, 30
fbounds = permutedims([extrema(ynoisy[.!isnan.(ynoisy)])...])
ymin, ymax = extrema(y)
fbounds[1] > ymin && (fbounds[1] = ymin)
fbounds[2] < ymax && (fbounds[2] = ymax)
sdev_prop = 0.05*diff(fbounds, dims=2)[:]
sdev_pos = [0.02abs(diff([extrema(x)...])[1])]
demean_ns = false
sampledc = true
K = transD_GP.GP.Mat32()
δ = 0.05
Random.seed!(13)
opt = transD_GP.OptionsNonstat(optlog10λ,
                        nmin = nmin,
                        nmax = nmax,
                        fbounds = fbounds,
                        δ = δ,
                        demean = demean_ns,
                        sampledc = sampledc,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        pnorm = pnorm,
                        K = K,
                        )
##

optlog10λ.history_mode = "a"
opt.history_mode = "a"
##
transD_GP.init_out_filenames(optlog10λ)
transD_GP.init_out_filenames(opt)
#bodge this so we don't have to run from the same directory as the old files
optlog10λ.fstar_filename = "testing_data/"*optlog10λ.fstar_filename
optlog10λ.x_ftrain_filename = "testing_data/"*optlog10λ.x_ftrain_filename
optlog10λ.costs_filename = "testing_data/"*optlog10λ.costs_filename
opt.fstar_filename = "testing_data/"*opt.fstar_filename
opt.x_ftrain_filename = "testing_data/"*opt.x_ftrain_filename
opt.costs_filename = "testing_data/"*opt.costs_filename

##
# correct values per chain, as stored in the test files
x_ftrain_ns = [[0.3606203901005928 0.40668939952022337 0.35695716468236716 0.2624710257822056 0.7309837184997988 0.09907569005826172 0.4384855047127546 0.4240562583253382 0.20723312563292223 0.0804216767745144 0.9006985134811043 0.3951757049312128 0.6269940818565283 0.058033000122784194 0.877004345217369 0.04005771809421885 0.6019455220987492 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
               -3.273539110218971 2.92765271708093 -2.8478724688082155 -0.6298259597050568 -1.4785412012542083 -1.8105167503237656 -4.2477591368662875 2.0208313246119136 -0.7930843122304994 -3.1330973202832118 2.1656551000180655 -1.3095036313557782 0.6874213494471846 -2.342087118177461 -4.146686443903445 -4.3536683926250435 -1.9897322978345238 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0],
              [0.41798164671939503 0.7476709354675064 0.3661786839239984 0.6298237044448398 0.11201264100758644 0.7285047725142729 0.46451120399414747 0.80641347881773 0.0960184161219579 0.7056667887095189 0.3582612012868851 0.6307055901719155 0.13157024627931138 0.6913116346109218 0.29197458639640855 0.10576318120807768 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
               2.9764886671528594 -1.5784857416346867 -0.85733122649682 -0.20058017471976688 -3.1734204772541665 0.38128726540862257 -3.780041779721524 3.4316019288954704 2.4754434998460173 -2.4391604603672983 0.6005511020333341 -3.695057348802962 3.913539538009535 2.1900600920728346 2.587739571040321 1.341864313226683 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0], 
              [0.37743600612320477 0.3985740427572615 0.30362739858964366 0.8809510167900401 0.4739648026251868 0.9196538218211477 0.775544477630246 0.4117556846060175 0.8344239816887886 0.7844413529690875 0.817627936701667 0.3984449091490326 0.3362651150450794 0.6037259556181453 0.4739347484294692 0.2003893167318971 0.6105120166255141 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
              -2.00988747151336 2.665474463799198 -0.7199830960773586 0.2146390981872257 2.461148229800755 0.583574636489491 -1.1969806846848172 -4.072215748620355 -3.489652641518525 2.9492819468178775 -4.396522341997929 -1.133203024189191 3.9634493307054983 -0.5926599413662617 -4.2296707991996785 -0.7095150807190489 1.2739513208331257 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0],
              [0.27869633855370896 0.9416066493267048 0.27940921149686404 0.5214656214322371 0.3858836080444921 0.18955378172250362 0.40307901488093106 0.12313470318106894 0.9418163392292713 0.12481948703853674 0.6768585318249375 0.19732460974737429 0.9723483772994661 0.4678583058972649 0.8426824556016561 0.34396040885327284 0.24488844225114917 0.5112911161047697 0.7745083466082023 0.09129437463075364 0.36892906140999937 0.5935645157732607 0.7488598975926185 0.5957267429123798 0.1630683809816259 0.0 0.0 0.0 0.0 0.0;
              -4.497476822329708 1.9121603990643552 2.4640232119706527 1.2616741887409304 -2.9423596767072278 -0.06857534330422055 2.3493265345553485 1.1757276273834814 3.090856955878926 -4.518770336744846 0.9472434204181379 -2.197223977042386 0.6488586458444761 2.6533975409629136 3.084231042792335 0.11185358358948072 -0.6780721519910976 -3.271241655225382 1.7491383780388192 1.6652422711460373 -0.7924652039627609 0.031137106251686697 3.34744423008942 -0.3308176518730255 -0.8258033014113981 0.0 0.0 0.0 0.0 0.0]]

n_ns = [4, 3, 6, 8]

x_ftrain_s = [[0.5026815805596105 0.5545227415572486 0.5146579819844821 0.37547140948856955 0.7890495747421975 0.18285276503217973 0.2516297666056329 0.9204359116840197 0.6799781601444546 0.8948396138231055 0.6229722184791631 0.1177948611059579 0.2310396863411785 0.9757606221577078 0.9200244653336747 0.3476422582739468 0.1884095848469977 0.28745982220232363 0.4861163927728539 0.3954356085167403 0.8464910848389889 0.2718982043089137 0.7727370785004052 0.8005979954704763 0.9774237218452777 0.9464807677622769 0.21934979776404562 0.005102141888382532 0.20925341847802229 0.4435978449911059; 
               -1.050559147483305 -0.8945196801924427 -0.9224293333195928 -0.9908148140709856 -0.9801577172088128 -0.9903442929041689 -1.0511384293647204 -1.059328707717616 -0.9591165125321414 -0.9947644414115122 -1.1835994371612621 -0.9578654658368201 -1.1364579111539308 -1.1244844474877818 -1.0186907511103582 -0.918044819019293 -0.8271972364142561 -0.8278438358663714 -1.1585506021823082 -1.1098064522373197 -1.0511904783374335 -1.077312065418988 -1.0463788590556076 -0.8563728031263591 -0.9929909748439182 -1.1965634457456695 -0.9520368953658574 -0.9671899444490628 -0.8206853753540002 -1.153507872014128],
              [0.8584134387357824 0.5464667199608733 0.3943588907839589 0.6847666793989239 0.20758593371584508 0.8460446095349564 0.12350050911292819 0.6046233505163399 0.8328630208149138 0.9526593337971843 0.007647097079578627 0.28410716566930166 0.17428369570156743 0.5119168650700942 0.3411122827959544 0.4069175792028222 0.7739240836283039 0.8080370255126386 0.23918364769075245 0.7961378311852106 0.24363167249638112 0.5436882883099681 0.037190642666712764 0.17542651911766355 0.13476666119609704 0.8646665158284379 0.3785179313182144 0.4216616345656411 0.8999185672516342 0.6023842252271996;
               -0.8209940512358518 -1.0798599334037402 -0.9352649178753314 -1.0202468736611556 -0.8541602836918148 -0.9203975522061871 -0.9156137352392366 -0.8580557108689477 -0.8703364533639802 -0.9713536758102499 -0.8441233801897189 -0.9984502140255286 -1.043869638677034 -1.197255300304413 -1.0267912484846842 -1.0690298677416974 -0.949114254302214 -0.8379685269119636 -1.0370417709011668 -0.914629154818752 -0.8863226803119202 -0.949216650228916 -1.1366003080978082 -0.8240851944606152 -0.8485441912965284 -1.0033651445970335 -1.0156379244892408 -1.0650564656877997 -1.1359443377242064 -1.1016135617759983], 
              [0.9605128272353466 0.38747513279645296 0.25191835762084436 0.23448969020267957 0.2725360886937186 0.22999848570255554 0.820227817588314 0.29182949628310906 0.340144244592684 0.5699779086957768 0.7873018783479561 0.5134512821655564 0.27118262493131645 0.11011724159867244 0.5295539583768022 0.1658120623973945 0.7116368653181858 0.6070172516202845 0.8645080557303385 0.05803511928493063 0.4812404701522854 0.4839386979645044 0.9006946924414921 0.6716038092376229 0.24410165819784937 0.4135248708111979 0.22806025837769936 0.26838028427348365 0.5576281711007943 0.676334452614073; 
               -1.096337414682399 -1.1792743083892825 -0.8641270356818355 -1.1426950568895016 -1.018806744129842 -1.1720674664137478 -0.81683895916483 -0.8962068699659175 -1.094516950356651 -0.9260612326616192 -1.1561078109543668 -1.0665676910528472 -0.9457120508446151 -1.0902981602517494 -1.0873626762665443 -1.045390496548077 -0.939014226047308 -0.9862832513473462 -1.1834335694270326 -0.9067413816103946 -0.9654805215474761 -1.0882408720437775 -1.1194426557464943 -1.0996219129845874 -1.0690965127432293 -0.8701522951515921 -0.9272630432793121 -1.0156292812242091 -0.9667039373733762 -0.8119865535470194], 
              [0.8156241624229966 0.44059562864761154 0.07157565635030176 0.20960634704047715 0.7405283291431856 0.033238442906150446 0.4623457448423742 0.4023764924977189 0.5902638619821537 0.39184175281314376 0.09476270926137245 0.42328732224788596 0.7041634503541843 0.5786709153831461 0.5628475368500884 0.7948437760956024 0.6910615021450255 0.8854924888027604 0.9084754812354597 0.9414963097210474 0.13929696003059985 0.12362240850218464 0.5324963446130311 0.5355023061561756 0.6122634387038886 0.9501580901245366 0.48458204159245094 0.9533890624956989 0.7889995211731183 0.5836750357218413;
               -0.8182522337085651 -1.0055213570763386 -1.0675537093092826 -0.9625072368228749 -0.9425486400402971 -0.8797404279363393 -1.0313965937732337 -1.1178111683981233 -0.8163497164903469 -1.038423348906223 -1.1830465282508673 -0.9054378597660917 -0.8467084151751388 -0.9988668266106078 -0.925461048210355 -0.9545108299723822 -0.9933631533886832 -0.9916676496038269 -1.039165147425735 -0.8617614433806986 -1.0053781609867825 -0.9958569196262486 -1.0614459751394383 -0.8177109218190248 -0.8668489574109816 -0.8207853651699865 -1.0944935575243915 -0.8041156327559184 -0.8808073749390615 -0.9602339184649771]]
n_s = [28, 25, 21, 23]

##
@testset "test re-initialising models from a file" begin
    models_stat = []
    models_nonstat = []
    @testset "initialise stationary model" begin
        for cidx = 1:4
            push!(models_stat, transD_GP.init(optlog10λ, cidx))
        end
        for cidx = 1:4
            @test models_stat[cidx].n == n_s[cidx]
            @test all(models_stat[cidx].xtrain[1:n_s[cidx]] .== x_ftrain_s[cidx][1,1:n_s[cidx]])
            @test all(models_stat[cidx].ftrain[1:n_s[cidx]] .== x_ftrain_s[cidx][2,1:n_s[cidx]])
        end
    end

    @testset "initialise nonstationary model" begin
        for cidx = 1:4
            push!(models_nonstat, transD_GP.init(opt, models_stat[cidx], cidx))
        end
        for cidx = 1:4
            @test models_nonstat[cidx].n == n_ns[cidx]
            @test all(models_nonstat[cidx].xtrain[1:n_ns[cidx]] .== x_ftrain_ns[cidx][1,1:n_ns[cidx]])
            @test all(models_nonstat[cidx].ftrain[1:n_ns[cidx]] .== x_ftrain_ns[cidx][2,1:n_ns[cidx]])
        end
    end

end
/*
RLFS = RIZ+linker+REZ
RIZ: GGGNGGGNGGG, each G cluster has at least 3 contiguous Gs, 1~2nt between clusters, at least three clusters
linker: 0~50nt
REZ: 100~2000nt, 40% G
*/

package main

import (
    "os"
    "io/ioutil"
    "fmt"
    "log"
    "net/http"
    "strings"
    "strconv"
    "regexp"
    "flag"
    "time"
)

var monitor = flag.Bool("m",false,"Monitor web requests")
var port = flag.Int("p",8686,"Server port")
var seqfile = flag.String("seqfile","","Sequence file")

const (
    version = "0.1.3"
    pageTop    = `<!DOCTYPE HTML><html><head>
<title>RLFS Finder</title>
<style>
.error{color:#FF0000;}
.riz{background-color:#FF0000;color:white}
.linker{background-color:#00FF00;color:white}
.rez{background-color:#5555FF;color:white}
.axis{color:#AAAAAA}
textarea{font-family:monospace}
</style>
</head>
<body>
<h3>RLFS Finder</h3>
<p><a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3258121/" target="_new">Reference</a>: Wongsurawat et al. Quantitative model of R-loop forming structures reveals a novel level of RNA-DNA interactome complexity. Nucleic Acids Res 40, e16 (2012)</p>
<p>Model: <b>RLFS = <span class="riz">RIZ</span>+<span class="linker">linker</span>+<span class="rez">REZ</span></b><p>`

    form       = `<form action="/" method="POST">
<label for="numbers">Sequence:</label><br />
<textarea name="sequence" rows="20" cols="100"></textarea>
<br/>
<input type="submit" value="Search">&nbsp;<input type="reset" value="Reset">
<br/>
</form>`

    pageBottom = `</body></html>`

    anError    = `<p class="error">%s</p>`
)

func main() {
    flag.Parse()
    if len(*seqfile) > 0{ //commandline mode, skipping web server
        dat, err := ioutil.ReadFile(*seqfile)
        if err != nil {
            panic(err)
        }
        seq := string(dat)
        segResults,_ := processSeqString(preProcessSeqString(seq))
        fmt.Println(len(segResults))
        os.Exit(0)
    }
    // Web server
    fmt.Println("RLFS Finder v"+version)
    if *port < 1 || *port > 65535 {
        log.Fatal("Invalid port ",*port)
    }
    var url string = "http://127.0.0.1:"+strconv.Itoa(*port)
    fmt.Println("Open link "+url+" in a web browser to access the web interface of this program")
    if err := openURL(url); err != nil {
        log.Println("Unable to open link in browser")
    }
    http.HandleFunc("/", homePage)
    if err := http.ListenAndServe(":"+strconv.Itoa(*port), nil); err != nil {
        log.Fatal("failed to start server", err)
    }
}

func homePage(writer http.ResponseWriter, request *http.Request) {
    err := request.ParseForm() // Must be called before writing response
    fmt.Fprint(writer, pageTop, form)
    if err != nil {
        fmt.Fprintf(writer, anError, err)
    } else {
        if seq, message, ok := processRequest(request); ok {
            preProcessedSeq := preProcessSeqString(seq)
            segResults,found := processSeqString(preProcessedSeq)
            if found {
                for i,v := range(segResults) {
                    fmt.Fprint(writer, "<h4>Match "+strconv.Itoa(i+1)+"</h4><pre>"+formatSeq(annotateSeq(v[0],v[1],v[2],v[3],v[4]),100)+"</pre><hr/>\n")
                }
            } else {
                fmt.Fprint(writer, "<h4>Pattern not found in sequence:</h4><pre>"+formatSeq(breakSeq(preProcessedSeq),100)+"</pre>\n")
            }
        } else if message != "" {
            fmt.Fprintf(writer, anError, message)
        }
    }
    fmt.Fprint(writer, pageBottom)
}

func processRequest(request *http.Request) (string,string,bool) {
    var seq string
    if *monitor {
        fmt.Println(time.Now(),request.RemoteAddr, request.UserAgent())
    }
    if slice, found := request.Form["sequence"]; found && len(slice) > 0 {
        seq = slice[0]
    }
    if len(seq) == 0 {
        return seq, "", false // no data first time form is shown
    }
    return seq, "", true
}

func preProcessSeqString(seq string) (string){
    seq = strings.ToLower(seq)
    seq = strings.Replace(seq," ","",-1)
    var newseq string
    for _,v := range(seq) {
        sv := string(v)
        if strings.Contains("acgt",sv) {
            newseq += sv
        }
    }
    return newseq
}

func processSeqString(seq string) ([][]string,bool) {
    rizResult,found := searchForRiz(seq)
    if !found {
        return nil,false
    }
    segSeqs := [][]string{}
    for _,riz := range(rizResult) {
        if len(riz[2]) < 100 {
            continue
        }
        seqLinker, seqRez, seqRest, found := searchForRez(riz[2])
        if !found {
            continue
        }
        segSeqs = append(segSeqs,[]string{riz[0],riz[1],seqLinker,seqRez,seqRest})
    }
    if len(segSeqs) == 0 {
        return nil,false
    } else {
        return segSeqs, true
    }
}

func searchForRiz(seq string) ([][]string,bool) {
    re := regexp.MustCompile("(g{3,}[act]{1,2}){2,}g{3,}")
    match := re.FindAllStringIndex(seq,-1)
    if match == nil {
        return nil,false
    }
    segSeqs := [][]string{}
    for _,m := range(match) {
        seqBeforeRiz := seq[:m[0]]
        seqRest := seq[m[1]:]
        var seqRiz string
        for _,v := range(seq[m[0]:m[1]]) {
            sv := string(v)
            seqRiz += sv
        }
        segSeqs = append(segSeqs,[]string{seqBeforeRiz,seqRiz,seqRest})
    }
    return segSeqs, true
}

func searchForRez(seq string) (string,string,string,bool) {
    seqLength := len(seq)
    for i := 0; i < 50; i++ {
        if (seqLength - i) < 100 {
            break
        }
        for j := 100; j <= 2000; j++ {
            if (seqLength - i - j) < 0 {
                break
            }
            if enoughG(seq[i:i+j]) {
                return seq[:i],seq[i:(i+j)],seq[(i+j):],true
            }
        }
    }
    return "","","",false
}

func enoughG(seq string) (bool) {
    return float64(strings.Count(seq,"g"))/float64(len(seq)) >= 0.4
}

func annotateSeq(seqBeforeRiz string,seqRiz string,seqLinker string,seqRez string,seqRest string) ([]string) {
    sliceSeq := []string{}
    var sv string
    for _,v := range(seqBeforeRiz) {
        sliceSeq = append(sliceSeq, string(v))
    }
    for _,v := range(seqRiz) {
        sv = string(v)
        if sv == "g" {
            sv = "G"
        }
        sliceSeq = append(sliceSeq, "<span class=\"riz\">"+sv+"</span>")
    }
    for _,v := range(seqLinker) {
        sliceSeq = append(sliceSeq, "<span class=\"linker\">"+string(v)+"</span>")
    }
    for _,v := range(seqRez) {
        sv = string(v)
        if sv == "g" {
            sv = "G"
        }
        sliceSeq = append(sliceSeq, "<span class=\"rez\">"+sv+"</span>")
    }
    for _,v := range(seqRest) {
        sliceSeq = append(sliceSeq, string(v))
    }
    return sliceSeq
}

func breakSeq(seq string) ([]string) {
    sliceSeq := []string{}
    for _,v := range(seq) {
        sliceSeq = append(sliceSeq,string(v))
    }
    return sliceSeq
}

func formatSeq(sliceSeq []string, basePerLine int) (string) {
    var outString string
    outString += "<span class=\"axis\">"
    for i := 1; i <= basePerLine; i++ {
        outString += strconv.Itoa(i % 10)
    }
    outString += "</span>"
    outString += "\n\n"
    for i := range(sliceSeq) {
        outString += sliceSeq[i]
        if (i + 1) % basePerLine == 0 {
            outString += "   <span class=\"axis\">"+strconv.Itoa(i+1)+"</span>\n"
        }
    }
    outString += "\n\n"
    outString += "<span class=\"axis\">"
    for i := 1; i <= basePerLine; i++ {
        outString += strconv.Itoa(i % 10)
    }
    outString += "</span>"

    return outString
}

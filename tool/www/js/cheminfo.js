function _classCallCheck(e, t) { if (!(e instanceof t)) throw new TypeError("Cannot call a class as a function") }
var _typeof = "function" == typeof Symbol && "symbol" == typeof Symbol.iterator ? function(e) { return typeof e } : function(e) { return e && "function" == typeof Symbol && e.constructor === Symbol && e !== Symbol.prototype ? "symbol" : typeof e },
    _createClass = function() {
        function e(e, t) { for (var i = 0; i < t.length; i++) { var r = t[i];
                r.enumerable = r.enumerable || !1, r.configurable = !0, "value" in r && (r.writable = !0), Object.defineProperty(e, r.key, r) } } return function(t, i, r) { return i && e(t.prototype, i), r && e(t, r), t } }(),
    ArrayHelper = function() {
        function e() { _classCallCheck(this, e) } return _createClass(e, null, [{ key: "clone", value: function(t) { var i = Array.isArray(t) ? [] : {}; for (var r in t) { var n = t[r]; "function" == typeof n.clone ? i[r] = n.clone() : i[r] = "object" === (void 0 === n ? "undefined" : _typeof(n)) ? e.clone(n) : n } return i } }, { key: "print", value: function(e) { if (0 == e.length) return ""; for (var t = "(", i = 0; i < e.length; i++) t += e[i].id ? e[i].id + ", " : e[i] + ", "; return (t = t.substring(0, t.length - 2)) + ")" } }, { key: "each", value: function(e, t) { for (var i = 0; i < e.length; i++) t(e[i]) } }, { key: "get", value: function(e, t, i) { for (var r = 0; r < e.length; r++)
                    if (e[r][t] == i) return e[r] } }, { key: "contains", value: function(e, t) { if (t.property || t.func) { if (t.func) { for (var i = 0; i < e.length; i++)
                            if (t.func(e[i])) return !0 } else
                        for (var r = 0; r < e.length; r++)
                            if (e[r][t.property] == t.value) return !0 } else
                    for (var n = 0; n < e.length; n++)
                        if (e[n] == t.value) return !0; return !1 } }, { key: "intersection", value: function(e, t) { for (var i = new Array, r = 0; r < e.length; r++)
                    for (var n = 0; n < t.length; n++) e[r] === t[n] && i.push(e[r]); return i } }, { key: "unique", value: function(e) { var t = {}; return e.filter(function(e) { return void 0 === t[e] && (t[e] = !0) }) } }, { key: "count", value: function(e, t) { for (var i = 0, r = 0; r < e.length; r++) e[r] === t && i++; return i } }, { key: "toggle", value: function(e, t) { for (var i = [], r = !1, n = 0; n < e.length; n++) e[n] !== t ? i.push(e[n]) : r = !0; return r || i.push(t), i } }, { key: "remove", value: function(e, t) { for (var i = [], r = 0; r < e.length; r++) e[r] !== t && i.push(e[r]); return i } }, { key: "removeAll", value: function(e, t) { return e.filter(function(e) { return -1 === t.indexOf(e) }) } }, { key: "merge", value: function(e, t) { for (var i = new Array(e.length + t.length), r = 0; r < e.length; r++) i[r] = e[r]; for (var n = 0; n < t.length; n++) i[e.length + n] = t[n]; return i } }, { key: "containsAll", value: function(e, t) { for (var i = 0, r = 0; r < e.length; r++)
                    for (var n = 0; n < t.length; n++) e[r] === t[n] && i++; return i === t.length } }, { key: "sortByAtomicNumberDesc", value: function(e) { var t = e.map(function(e, t) { return { index: t, value: e.atomicNumber.split(".").map(Number) } }); return t.sort(function(e, t) { for (var i = Math.min(t.value.length, e.value.length), r = 0; r < i && t.value[r] === e.value[r];) r++; return r === i ? t.value.length - e.value.length : t.value[r] - e.value[r] }), t.map(function(t) { return e[t.index] }) } }]), e }(),
    Atom = function() {
        function e(t) { var i = arguments.length > 1 && void 0 !== arguments[1] ? arguments[1] : "-";
            _classCallCheck(this, e), this.element = t, this.explicit = !1, this.ringbonds = new Array, this.rings = new Array, this.bondType = i, this.isTerminal = !1, this.isBridge = !1, this.isBridgeNode = !1, this.originalRings = new Array, this.bridgedRing = null, this.anchoredRings = new Array, this.bracket = null, this.chiral = 0, this.order = {} } return _createClass(e, [{ key: "addAnchoredRing", value: function(e) { ArrayHelper.contains(this.anchoredRings, { value: e }) || this.anchoredRings.push(e) } }, { key: "getRingbondCount", value: function() { return this.ringbonds.length } }, { key: "canRotate", value: function() { return "-" === this.bondType && 0 == this.rings.length } }, { key: "hasRingbonds", value: function() { return this.ringbonds.length > 0 } }, { key: "getMaxRingbond", value: function() { for (var e = 0, t = 0; t < this.ringbonds.length; t++) this.ringbonds[t].id > e && (e = this.ringbonds[t].id); return e } }, { key: "hasRing", value: function(e) { for (var t = 0; t < this.rings; t++)
                    if (e === this.rings[t]) return !0;
                return !1 } }, { key: "backupRings", value: function() { this.originalRings = []; for (var e = 0; e < this.rings.length; e++) this.originalRings.push(this.rings[e]) } }, { key: "restoreRings", value: function() { this.rings = []; for (var e = 0; e < this.originalRings.length; e++) this.rings.push(this.originalRings[e]) } }, { key: "haveCommonRingbond", value: function(e, t) { for (var i = 0; i < e.ringbonds.length; i++)
                    for (var r = 0; r < t.ringbonds.length; r++)
                        if (e.ringbonds[i].id == t.ringbonds[r].id) return !0;
                return !1 } }, { key: "maxCommonRingbond", value: function(e, t) { for (var i = 0, r = 0, n = 0, s = 0; s < e.ringbonds.length; s++) { e.ringbonds[s].id > r && (r = e.ringbonds[s].id); for (var o = 0; o < t.ringbonds.length; o++) t.ringbonds[o].id > n ? n = t.ringbonds[o].id : r == n && (i = r) } return i } }, { key: "getOrder", value: function(e) { return this.order[e] } }, { key: "setOrder", value: function(e, t) { this.order[e] = t } }, { key: "getAtomicNumber", value: function() { return e.atomicNumbers[this.element] } }], [{ key: "sortByAtomicNumber", value: function(t, i) { for (var r = new Array(t.length), n = 0; n < t.length; n++) { var s = i[t[n]],
                        o = e.atomicNumbers[s.value.element];
                    r[n] = { atomicNumber: o.toString(), vertexId: s.id } } return ArrayHelper.sortByAtomicNumberDesc(r) } }, { key: "hasDuplicateAtomicNumbers", value: function(e) { for (var t = {}, i = 0; i < e.length; i++) { var r = e[i]; if (void 0 !== t[r.atomicNumber]) return !0;
                    t[r.atomicNumber] = !0 } return !1 } }, { key: "getDuplicateAtomicNumbers", value: function(e) { for (var t = {}, i = [], r = 0; r < e.length; r++) { var n = e[r];
                    void 0 === t[n.atomicNumber] && (t[n.atomicNumber] = []), t[n.atomicNumber].push(r) } for (var s in t) { var o = t[s];
                    o.length > 1 && i.push(o) } return i } }]), e }();
Atom.atomicNumbers = { H: 1, He: 2, Li: 3, Be: 4, B: 5, b: 5, C: 6, c: 6, N: 7, n: 7, O: 8, o: 8, F: 9, Ne: 10, Na: 11, Mg: 12, Al: 13, Si: 14, P: 15, p: 15, S: 16, s: 16, Cl: 17, Ar: 18, K: 19, Ca: 20, Sc: 21, Ti: 22, V: 23, Cr: 24, Mn: 25, Fe: 26, Co: 27, Ni: 28, Cu: 29, Zn: 30, Ga: 31, Ge: 32, As: 33, Se: 34, Br: 35, Kr: 36, Rb: 37, Sr: 38, Y: 39, Zr: 40, Nb: 41, Mo: 42, Tc: 43, Ru: 44, Rh: 45, Pd: 46, Ag: 47, Cd: 48, In: 49, Sn: 50, Sb: 51, Te: 52, I: 53, Xe: 54, Cs: 55, Ba: 56, La: 57, Ce: 58, Pr: 59, Nd: 60, Pm: 61, Sm: 62, Eu: 63, Gd: 64, Tb: 65, Dy: 66, Ho: 67, Er: 68, Tm: 69, Yb: 70, Lu: 71, Hf: 72, Ta: 73, W: 74, Re: 75, Os: 76, Ir: 77, Pt: 78, Au: 79, Hg: 80, Tl: 81, Pb: 82, Bi: 83, Po: 84, At: 85, Rn: 86, Fr: 87, Ra: 88, Ac: 89, Th: 90, Pa: 91, U: 92, Np: 93, Pu: 94, Am: 95, Cm: 96, Bk: 97, Cf: 98, Es: 99, Fm: 100, Md: 101, No: 102, Lr: 103, Rf: 104, Db: 105, Sg: 106, Bh: 107, Hs: 108, Mt: 109, Ds: 110, Rg: 111, Cn: 112, Uut: 113, Uuq: 114, Uup: 115, Uuh: 116, Uus: 117, Uuo: 118 };
var CanvasWrapper = function() {
        function e(t, i, r, n) { _classCallCheck(this, e), "string" == typeof t || t instanceof String ? this.canvas = document.getElementById(t) : this.canvas = t, this.ctx = this.canvas.getContext("2d"), this.colors = i, this.bondLength = r, this.bondSpacing = n, this.drawingWidth = 0, this.drawingHeight = 0, this.offsetX = 0, this.offsetY = 0, this.clear() } return _createClass(e, [{ key: "setTheme", value: function(e) { this.colors = e } }, { key: "scale", value: function(e) { for (var t = { x: -Number.MAX_VALUE, y: -Number.MAX_VALUE }, i = { x: Number.MAX_VALUE, y: Number.MAX_VALUE }, r = 0; r < e.length; r++) { var n = e[r].position;
                    t.x < n.x && (t.x = n.x), t.y < n.y && (t.y = n.y), i.x > n.x && (i.x = n.x), i.y > n.y && (i.y = n.y) }
                t.x += 20, t.y += 20, i.x -= 20, i.y -= 20, this.drawingWidth = t.x - i.x, this.drawingHeight = t.y - i.y; var s = this.canvas.offsetWidth / this.drawingWidth,
                    o = this.canvas.offsetHeight / this.drawingHeight,
                    a = s < o ? s : o;
                this.ctx.scale(a, a), this.offsetX = -i.x, this.offsetY = -i.y, s < o ? this.offsetY += this.canvas.offsetHeight / (2 * a) - this.drawingHeight / 2 : this.offsetX += this.canvas.offsetWidth / (2 * a) - this.drawingWidth / 2 } }, { key: "reset", value: function() { this.ctx.setTransform(1, 0, 0, 1, 0, 0) } }, { key: "getColor", value: function(e) { return e = e.toUpperCase(), e in this.colors ? this.colors[e] : this.colors.C } }, { key: "drawCircle", value: function(e, t, i, r) { var n = !(arguments.length > 4 && void 0 !== arguments[4]) || arguments[4],
                    s = arguments.length > 5 && void 0 !== arguments[5] && arguments[5],
                    o = arguments.length > 6 && void 0 !== arguments[6] ? arguments[6] : ""; if (!isNaN(e) && !isNaN(t)) { var a = this.ctx,
                        l = this.offsetX,
                        h = this.offsetY;
                    a.save(), a.lineWidth = 1.5, a.beginPath(), a.arc(e + l, t + h, i, 0, 2 * Math.PI, !0), a.closePath(), s ? (n ? (a.fillStyle = "#f00", a.fill()) : (a.strokeStyle = "#f00", a.stroke()), this.drawDebugText(e, t, o)) : n ? (a.fillStyle = r, a.fill()) : (a.strokeStyle = r, a.stroke()), a.restore() } } }, { key: "drawLine", value: function(e) { if (!(isNaN(e.from.x) || isNaN(e.from.y) || isNaN(e.to.x) || isNaN(e.to.y))) { var t = this.ctx,
                        i = this.offsetX,
                        r = this.offsetY,
                        n = e.clone().shorten(8),
                        s = n.getLeftVector().clone(),
                        o = n.getRightVector().clone();
                    s.x += i, s.y += r, o.x += i, o.y += r, s = e.getLeftVector().clone(), o = e.getRightVector().clone(), s.x += i, s.y += r, o.x += i, o.y += r, t.save(), t.beginPath(), t.moveTo(s.x, s.y), t.lineTo(o.x, o.y), t.lineCap = "round", t.lineWidth = 1.5; var a = this.ctx.createLinearGradient(s.x, s.y, o.x, o.y);
                    a.addColorStop(.4, this.getColor(e.getLeftElement()) || this.getColor("C")), a.addColorStop(.6, this.getColor(e.getRightElement()) || this.getColor("C")), t.strokeStyle = a, t.stroke(), t.restore() } } }, { key: "drawWedge", value: function(e) { var t = arguments.length > 1 && void 0 !== arguments[1] ? arguments[1] : 3; if (!(isNaN(e.from.x) || isNaN(e.from.y) || isNaN(e.to.x) || isNaN(e.to.y))) { var i = this.ctx,
                        r = this.offsetX,
                        n = this.offsetY,
                        s = e.clone().shorten(8),
                        o = s.getLeftVector().clone(),
                        a = s.getRightVector().clone();
                    o.x += r, o.y += n, a.x += r, a.y += n, o = e.getLeftVector().clone(), a = e.getRightVector().clone(), o.x += r, o.y += n, a.x += r, a.y += n, i.save(); var l = Vector2.normals(o, a);
                    l[0].normalize(), l[1].normalize(); var h = e.getRightChiral(),
                        u = o,
                        g = a;
                    h && (u = a, g = o); var c = Vector2.add(u, Vector2.multiply(l[0], .75)),
                        v = Vector2.add(g, Vector2.multiply(l[0], t)),
                        d = Vector2.add(g, Vector2.multiply(l[1], t)),
                        f = Vector2.add(u, Vector2.multiply(l[1], .75));
                    i.beginPath(), i.moveTo(c.x, c.y), i.lineTo(v.x, v.y), i.lineTo(d.x, d.y), i.lineTo(f.x, f.y); var p = this.ctx.createRadialGradient(a.x, a.y, this.bondLength, a.x, a.y, 0);
                    p.addColorStop(.4, this.getColor(e.getLeftElement()) || this.getColor("C")), p.addColorStop(.6, this.getColor(e.getRightElement()) || this.getColor("C")), i.fillStyle = p, i.fill(), i.restore() } } }, { key: "drawDashedWedge", value: function(e) { var t = arguments.length > 1 && void 0 !== arguments[1] ? arguments[1] : 6; if (!(isNaN(e.from.x) || isNaN(e.from.y) || isNaN(e.to.x) || isNaN(e.to.y))) { var i = this.ctx,
                        r = this.offsetX,
                        n = this.offsetY,
                        s = e.getLeftVector().clone(),
                        o = e.getRightVector().clone();
                    s.x += r, s.y += n, o.x += r, o.y += n, i.save(); var a = Vector2.normals(s, o);
                    a[0].normalize(), a[1].normalize(); var l = e.getRightChiral(),
                        h = void 0,
                        u = void 0,
                        g = void 0,
                        c = void 0,
                        v = e.clone();
                    l ? (h = o, u = s, v.shortenRight(3), g = v.getRightVector().clone(), c = v.getLeftVector().clone()) : (h = s, u = o, v.shortenLeft(3), g = v.getLeftVector().clone(), c = v.getRightVector().clone()), g.x += r, g.y += n, c.x += r, c.y += n; var d = Vector2.add(h, Vector2.multiply(a[0], .75)),
                        f = Vector2.add(u, Vector2.multiply(a[0], t / 2)),
                        p = Vector2.add(u, Vector2.multiply(a[1], t / 2)),
                        y = Vector2.add(h, Vector2.multiply(a[1], .75));
                    i.beginPath(), i.moveTo(d.x, d.y), i.lineTo(f.x, f.y), i.lineTo(p.x, p.y), i.lineTo(y.x, y.y); var m = this.ctx.createRadialGradient(o.x, o.y, this.bondLength, o.x, o.y, 0);
                    m.addColorStop(.4, this.getColor(e.getLeftElement()) || this.getColor("C")), m.addColorStop(.6, this.getColor(e.getRightElement()) || this.getColor("C")), i.fillStyle = m, i.fill(), i.globalCompositeOperation = "destination-out", i.beginPath(), i.moveTo(g.x, g.y), i.lineTo(c.x, c.y), i.lineCap = "butt", i.lineWidth = t, i.setLineDash([1, 1]), i.strokeStyle = this.getColor("BACKGROUND"), i.stroke(), i.globalCompositeOperation = "source-over", i.restore() } } }, { key: "drawDebugText", value: function(e, t, i) { var r = this.ctx;
                r.save(), r.font = "5px Droid Sans, sans-serif", r.textAlign = "start", r.textBaseline = "top", r.fillStyle = "#ff0000", r.fillText(i, e + this.offsetX, t + this.offsetY), r.restore() } }, { key: "drawBall", value: function(e, t, i) { var r = this.ctx;
                r.save(), r.beginPath(), r.arc(e + this.offsetX, t + this.offsetY, this.bondLength / 4.5, 0, 2 * Math.PI, !1), r.fillStyle = this.getColor(i), r.fill() } }, { key: "drawText", value: function(e, t, i, r, n, s, o) { if (!isNaN(e) && !isNaN(t)) { var a = this.ctx,
                        l = this.offsetX,
                        h = this.offsetY;
                    a.save(); var u = "10px Droid Sans, sans-serif",
                        g = "6px Droid Sans, sans-serif";
                    a.textAlign = "start", a.textBaseline = "top"; var c = "+",
                        v = 0;
                    o && (2 === o ? c = "2+" : -1 === o ? c = "-" : -2 === o && (c = "2-"), a.font = g, v = a.measureText(c).width), a.font = u, a.fillStyle = this.getColor(i); var d = a.measureText(i);
                    d.totalWidth = d.width + v, d.height = parseInt(u, 10); var f = d.totalWidth > d.height ? d.totalWidth : d.height;
                    f /= 2, a.globalCompositeOperation = "destination-out", a.beginPath(), a.arc(e + l, t + h + d.height / 20, f + 1, 0, 2 * Math.PI, !0), a.closePath(), a.fill(), a.globalCompositeOperation = "source-over", a.fillStyle = this.getColor(i), a.fillText(i, e - d.totalWidth / 2 + l, t - d.height / 2 + h), o && (a.font = g, a.fillText(c, e - d.totalWidth / 2 + d.width + l, t - d.height / 2 + h)), a.font = u; var p = a.measureText("H"); if (p.height = parseInt(u, 10), 1 === r) { var y = e - d.totalWidth / 2 + l,
                            m = t - d.height / 2 + h; "left" === n ? y -= d.totalWidth : "right" === n ? y += d.totalWidth : "up" === n && s ? y += d.totalWidth : "down" === n && s ? y += d.totalWidth : "up" !== n || s ? "down" !== n || s || (m += d.height) : m -= d.height, a.fillText("H", y, m) } else if (r > 1) { var b = e - d.totalWidth / 2 + l,
                            A = t - d.height / 2 + h;
                        a.font = g; var x = a.measureText(r);
                        x.height = parseInt(g, 10), "left" === n ? b -= p.width + x.width : "right" === n ? b += d.totalWidth : "up" === n && s ? b += d.totalWidth : "down" === n && s ? b += d.totalWidth : "up" !== n || s ? "down" !== n || s || (A += d.height) : A -= d.height, a.font = u, a.fillText("H", b, A), a.font = g, a.fillText(r, b + p.width, A + p.height / 2) }
                    a.restore() } } }, { key: "drawDebugPoint", value: function(e, t) { var i = arguments.length > 2 && void 0 !== arguments[2] ? arguments[2] : "",
                    r = arguments.length > 3 && void 0 !== arguments[3] ? arguments[3] : "#f00";
                this.drawCircle(e, t, 2, r, !0, !0, i) } }, { key: "drawAromaticityRing", value: function(e) { var t = this.ctx,
                    i = MathHelper.polyCircumradius(this.bondLength, e.getSize());
                t.save(), t.strokeStyle = this.getColor("C"), t.lineWidth = 1.5, t.beginPath(), t.arc(e.center.x + this.offsetX, e.center.y + this.offsetY, i - this.bondLength / 3 - this.bondSpacing, 0, 2 * Math.PI, !0), t.closePath(), t.stroke(), t.restore() } }, { key: "clear", value: function() { this.ctx.clearRect(0, 0, this.canvas.offsetWidth, this.canvas.offsetHeight) } }]), e }(),
    Edge = function() {
        function e(t, i, r) { _classCallCheck(this, e), this.id = null, this.sourceId = t, this.targetId = i, this.weight = r, this.bondType = "-", this.isInAromaticRing = !1, this.center = !1, this.chiral = "" } return _createClass(e, [{ key: "getBondCount", value: function() { return e.bonds[this.bondType] } }], [{ key: "bonds", get: function() { return { "-": 1, "/": 1, "\\": 1, "=": 2, "#": 3, $: 4 } } }]), e }(),
    Line = function() {
        function e() { var t = arguments.length > 0 && void 0 !== arguments[0] ? arguments[0] : new Vector2(0, 0),
                i = arguments.length > 1 && void 0 !== arguments[1] ? arguments[1] : new Vector(0, 0),
                r = arguments.length > 2 && void 0 !== arguments[2] ? arguments[2] : null,
                n = arguments.length > 3 && void 0 !== arguments[3] ? arguments[3] : null,
                s = arguments.length > 4 && void 0 !== arguments[4] && arguments[4],
                o = arguments.length > 5 && void 0 !== arguments[5] && arguments[5];
            _classCallCheck(this, e), this.from = t, this.to = i, this.elementFrom = r, this.elementTo = n, this.chiralFrom = s, this.chiralTo = o } return _createClass(e, [{ key: "clone", value: function() { return new e(this.from.clone(), this.to.clone(), this.elementFrom, this.elementTo) } }, { key: "getLength", value: function() { return Math.sqrt(Math.pow(this.to.x - this.from.x, 2) + Math.pow(this.to.y - this.from.y, 2)) } }, { key: "getAngle", value: function() { return Vector2.subtract(this.getRightVector(), this.getLeftVector()).angle() } }, { key: "getRightVector", value: function() { return this.from.x < this.to.x ? this.to : this.from } }, { key: "getLeftVector", value: function() { return this.from.x < this.to.x ? this.from : this.to } }, { key: "getRightElement", value: function() { return this.from.x < this.to.x ? this.elementTo : this.elementFrom } }, { key: "getLeftElement", value: function() { return this.from.x < this.to.x ? this.elementFrom : this.elementTo } }, { key: "getRightChiral", value: function() { return this.from.x < this.to.x ? this.chiralTo : this.chiralFrom } }, { key: "getLeftChiral", value: function() { return this.from.x < this.to.x ? this.chiralFrom : this.chiralTo } }, { key: "setRightVector", value: function(e, t) { return this.from.x < this.to.x ? this.to.set(e, t) : this.from.set(e, t), this } }, { key: "setLeftVector", value: function(e, t) { return this.from.x < this.to.x ? this.from.set(e, t) : this.to.set(e, t), this } }, { key: "rotateToXAxis", value: function() { var e = this.getLeftVector(); return this.setRightVector(e.x + this.getLength(), e.y), this } }, { key: "rotate", value: function(e) { var t = this.getLeftVector(),
                    i = this.getRightVector(),
                    r = Math.cos(e) * (i.x - t.x) - Math.sin(e) * (i.y - t.y) + t.x,
                    n = Math.sin(e) * (i.x - t.x) - Math.cos(e) * (i.y - t.y) + t.y; return this.setRightVector(r, n), this } }, { key: "shortenFrom", value: function(e) { var t = Vector2.subtract(this.to, this.from); return t.normalize(), t.multiply(e), this.from.add(t), this } }, { key: "shortenTo", value: function(e) { var t = Vector2.subtract(this.from, this.to); return t.normalize(), t.multiply(e), this.to.add(t), this } }, { key: "shortenRight", value: function(e) { return this.from.x < this.to.x ? this.shortenTo(e) : this.shortenFrom(e), this } }, { key: "shortenLeft", value: function(e) { return this.from.x < this.to.x ? this.shortenFrom(e) : this.shortenTo(e), this } }, { key: "shorten", value: function(e) { var t = Vector2.subtract(this.from, this.to); return t.normalize(), t.multiply(e / 2), this.to.add(t), this.from.subtract(t), this } }, { key: "getNormals", value: function() { var e = Vector2.subtract(from, to); return [new Vector2(-e.y, e.x), new Vector2(e.y, -e.x)] } }]), e }(),
    MathHelper = function() {
        function e() { _classCallCheck(this, e) } return _createClass(e, null, [{ key: "round", value: function(e, t) { return t = t || 1, Number(Math.round(e + "e" + t) + "e-" + t) } }, { key: "meanAngle", value: function(e) { for (var t = 0, i = 0, r = 0; r < e.length; r++) t += Math.sin(e[r]), i += Math.cos(e[r]); return Math.atan2(t / e.length, i / e.length) } }, { key: "innerAngle", value: function(t) { return e.toRad(180 * (t - 2) / t) } }, { key: "polyCircumradius", value: function(e, t) { return e / (2 * Math.sin(Math.PI / t)) } }, { key: "apothem", value: function(e, t) { return e * Math.cos(Math.PI / t) } }, { key: "apothemFromSideLength", value: function(t, i) { var r = e.polyCircumradius(t, i); return e.apothem(r, i) } }, { key: "centralAngle", value: function(t) { return e.toRad(360 / t) } }, { key: "toDeg", value: function(e) { return 180 * e / Math.PI } }, { key: "toRad", value: function(e) { return e * Math.PI / 180 } }]), e }(),
    Pair = function() {
        function e(t, i) { _classCallCheck(this, e), this.first = t, this.second = i } return _createClass(e, [{ key: "getHash", value: function() { return .5 * (this.first + this.second) * (this.first + this.second + 1) + (this.first + this.second) } }, { key: "contains", value: function(e) { return this.first === e || this.second === e } }], [{ key: "createUniquePairs", value: function(t) { for (var i = [], r = 0; r < t.length; r++)
                    for (var n = t[r], s = r + 1; s < t.length; s++) { var o = t[s];
                        i.push(new e(n, o)) }
                return i } }]), e }(),
    Ring = function() {
        function e(t, i, r) { _classCallCheck(this, e), this.id = null, this.ringbond = t, this.sourceId = i, this.targetId = r, this.members = new Array, this.edges = new Array, this.insiders = new Array, this.neighbours = new Array, this.positioned = !1, this.center = new Vector2, this.rings = new Array, this.isBridged = !1, this.template = null, this.isSpiro = !1, this.isFused = !1, this.centralAngle = 0, this.canFlip = !0 } return _createClass(e, [{ key: "clone", value: function() { var t = new e(this.ringbond, this.sourceId, this.target); return t.id = this.id, t.members = ArrayHelper.clone(this.members), t.insiders = ArrayHelper.clone(this.insiders), t.neighbours = ArrayHelper.clone(this.neighbours), t.positioned = this.positioned, t.center = this.center.clone(), t.rings = ArrayHelper.clone(this.rings), t.isBridged = this.isBridged, t.template = this.template, t.isSpiro = this.isSpiro, t.isFused = this.isFused, t.centralAngle = this.centralAngle, t.canFlip = this.canFlip, t } }, { key: "allowsFlip", value: function() { return this.canFlip && this.members.length > 4 } }, { key: "setFlipped", value: function() { this.members.length < 8 && (this.canFlip = !1) } }, { key: "getSize", value: function() { return this.members.length } }, { key: "getPolygon", value: function(e) { for (var t = [], i = 0; i < this.members.length; i++) t.push(e[this.members[i]].position); return t } }, { key: "getAngle", value: function() { return Math.PI - this.centralAngle } }, { key: "eachMember", value: function(e, t, i, r) { i = i || 0 === i ? i : this.members[0]; for (var n = i, s = 0; null != n && s < 100;) { var o = n; if (t(o), n = e[n].getNextInRing(e, this.id, r), r = o, n == i && (n = null), 99 == s) throw console.log("Smiles-drawer was not able to loop over the members of this ring.", this), "Smiles-drawer was not able to loop over the members of this ring.";
                    s++ } } }, { key: "getOrderedNeighbours", value: function(e) { for (var t = [], i = 0; i < this.neighbours.length; i++) { var r = RingConnection.getVertices(e, this.id, this.neighbours[i]);
                    t.push({ n: r.length, neighbour: this.neighbours[i] }) } return t.sort(function(e, t) { return t.n - e.n }), t } }, { key: "isAromatic", value: function(e) { for (var t = 0; t < this.members.length; t++) { var i = e[this.members[t]].value.element.charAt(0); if (i === i.toUpperCase()) return !1 } return !0 } }, { key: "isBenzeneLike", value: function(e) { var t = this.getDoubleBondCount(e),
                    i = this.members.length; return 3 === t && 6 === i || 2 === t && 5 === i } }, { key: "getDoubleBondCount", value: function(e) { for (var t = 0, i = 0; i < this.members.length; i++) { var r = e[this.members[i]].value; "=" !== r.bondType && "=" !== r.branchBond || t++ } return t } }, { key: "contains", value: function(e) { for (var t = 0; t < this.members.length; t++)
                    if (this.members[t] == e) return !0;
                return !1 } }, { key: "thisOrNeighboursContain", value: function(t, i) { for (var r = 0; r < this.neighbours.length; r++)
                    if (e.getRing(t, this.neighbours[r]).contains(i)) return !0;
                return !!this.contains(i) } }, { key: "hasSource", value: function() { return !(void 0 === this.sourceId || null === this.sourceId) } }, { key: "hasTarget", value: function() { return !(void 0 === this.targetId || null === this.targetId) } }, { key: "hasSourceAndTarget", value: function() { return this.hasSource() && this.hasTarget() } }], [{ key: "getRing", value: function(e, t) { for (var i = 0; i < e.length; i++)
                    if (e[i].id == t) return e[i] } }]), e }(),
    RingConnection = function() {
        function e(t, i) { _classCallCheck(this, e), this.id = null, this.rings = new Pair(t.id, i.id), this.vertices = []; for (var r = 0; r < t.members.length; r++)
                for (var n = t.members[r], s = 0; s < i.members.length; s++) { var o = i.members[s];
                    n === o && this.addVertex(n) } } return _createClass(e, [{ key: "addVertex", value: function(e) { ArrayHelper.contains(this.vertices, { value: e }) || this.vertices.push(e) } }, { key: "isBridge", value: function(e) { if (this.vertices.length > 2) return !0; for (var t = 0; t < this.vertices.length; t++) { if (e[this.vertices[t]].value.rings.length > 2) return !0 } return !1 } }, { key: "updateOther", value: function(e, t) { this.rings.first === t ? this.rings.second = e : this.rings.first = e } }], [{ key: "isBridge", value: function(e, t, i, r) { for (var n = null, s = 0; s < e.length; s++) { n = e[s]; var o = n.rings; if (o.first === i && o.second === r || o.first === r && o.second === i) return n.isBridge(t) } return !1 } }, { key: "getNeighbours", value: function(e, t) { for (var i = [], r = 0; r < e.length; r++) { var n = e[r].rings;
                    n.first === t ? i.push(n.second) : n.second === t && i.push(n.first) } return i } }, { key: "getVertices", value: function(e, t, i) { for (var r = 0; r < e.length; r++) { var n = e[r]; if (n.rings.first == t && n.rings.second == i || n.rings.first == i && n.rings.second == t) return n.vertices } } }]), e }(),
    SMILESPARSER = function() {
        function e(t, i, r, n) { this.message = t, this.expected = i, this.found = r, this.location = n, this.name = "SyntaxError", "function" == typeof Error.captureStackTrace && Error.captureStackTrace(this, e) }

        function t(t) {
            function i(e) { var i, r, n = it[e]; if (n) return n; for (i = e - 1; !it[i];) i--; for (n = it[i], n = { line: n.line, column: n.column, seenCR: n.seenCR }; i < e;) r = t.charAt(i), "\n" === r ? (n.seenCR || n.line++, n.column = 1, n.seenCR = !1) : "\r" === r || "\u2028" === r || "\u2029" === r ? (n.line++, n.column = 1, n.seenCR = !0) : (n.column++, n.seenCR = !1), i++; return it[e] = n, n }

            function r(e, t) { var r = i(e),
                    n = i(t); return { start: { offset: e, line: r.line, column: r.column }, end: { offset: t, line: n.line, column: n.column } } }

            function n(e) { et < rt || (et > rt && (rt = et, nt = []), nt.push(e)) }

            function s(t, i, r, n) { return null !== i && function(e) { var t = 1; for (e.sort(function(e, t) { return e.description < t.description ? -1 : e.description > t.description ? 1 : 0 }); t < e.length;) e[t - 1] === e[t] ? e.splice(t, 1) : t++ }(i), new e(null !== t ? t : function(e, t) { var i, r, n, s = new Array(e.length); for (n = 0; n < e.length; n++) s[n] = e[n].description; return i = e.length > 1 ? s.slice(0, -1).join(", ") + " or " + s[e.length - 1] : s[0], r = t ? '"' + function(e) {
                        function t(e) { return e.charCodeAt(0).toString(16).toUpperCase() } return e.replace(/\\/g, "\\\\").replace(/"/g, '\\"').replace(/\x08/g, "\\b").replace(/\t/g, "\\t").replace(/\n/g, "\\n").replace(/\f/g, "\\f").replace(/\r/g, "\\r").replace(/[\x00-\x07\x0B\x0E\x0F]/g, function(e) { return "\\x0" + t(e) }).replace(/[\x10-\x1F\x80-\xFF]/g, function(e) { return "\\x" + t(e) }).replace(/[\u0100-\u0FFF]/g, function(e) { return "\\u0" + t(e) }).replace(/[\u1000-\uFFFF]/g, function(e) { return "\\u" + t(e) }) }(t) + '"' : "end of input", "Expected " + i + " but " + r + " found." }(i, r), i, r, n) }

            function o() { var e, t, i, r, n, s, u, g, c, v; if (e = et, t = et, (i = l()) !== w) { for (r = [], n = a(); n !== w;) r.push(n), n = a(); if (r !== w) { for (n = [], s = et, u = h(), u === w && (u = null), u !== w ? (g = f(), g !== w ? (u = [u, g], s = u) : (et = s, s = w)) : (et = s, s = w); s !== w;) n.push(s), s = et, u = h(), u === w && (u = null), u !== w ? (g = f(), g !== w ? (u = [u, g], s = u) : (et = s, s = w)) : (et = s, s = w); if (n !== w) { for (s = [], u = a(); u !== w;) s.push(u), u = a(); if (s !== w)
                                if (u = h(), u === w && (u = null), u !== w)
                                    if (g = o(), g === w && (g = null), g !== w) { for (c = [], v = a(); v !== w;) c.push(v), v = a();
                                        c !== w ? (i = [i, r, n, s, u, g, c], t = i) : (et = t, t = w) } else et = t, t = w;
                            else et = t, t = w;
                            else et = t, t = w } else et = t, t = w } else et = t, t = w } else et = t, t = w; return t !== w && (tt = e, t = N(t)), e = t }

            function a() { var e, i, r, s, a, l; return e = et, i = et, 40 === t.charCodeAt(et) ? (r = L, et++) : (r = w, 0 === st && n(B)), r !== w ? (s = h(), s === w && (s = null), s !== w ? (a = o(), a !== w ? (41 === t.charCodeAt(et) ? (l = T, et++) : (l = w, 0 === st && n(I)), l !== w ? (r = [r, s, a, l], i = r) : (et = i, i = w)) : (et = i, i = w)) : (et = i, i = w)) : (et = i, i = w), i !== w && (tt = e, i = H(i)), e = i }

            function l() { var e, t; return e = et, t = g(), t === w && (t = c()) === w && (t = u()) === w && (t = v()), t !== w && (tt = e, t = M(t)), e = t }

            function h() { var e, i; return e = et, P.test(t.charAt(et)) ? (i = t.charAt(et), et++) : (i = w, 0 === st && n(O)), i !== w && (tt = e, i = W(i)), e = i }

            function u() { var e, i, r, s, o, a, l, h, u, g; return e = et, i = et, 91 === t.charCodeAt(et) ? (r = F, et++) : (r = w, 0 === st && n(E)), r !== w ? (s = k(), s === w && (s = null), s !== w ? (t.substr(et, 2) === D ? (o = D, et += 2) : (o = w, 0 === st && n(z)), o === w && (t.substr(et, 2) === _ ? (o = _, et += 2) : (o = w, 0 === st && n(q)), o === w && (o = c()) === w && (o = d()) === w && (o = v())), o !== w ? (a = p(), a === w && (a = null), a !== w ? (l = A(), l === w && (l = null), l !== w ? (h = y(), h === w && (h = null), h !== w ? (u = x(), u === w && (u = null), u !== w ? (93 === t.charCodeAt(et) ? (g = U, et++) : (g = w, 0 === st && n(j)), g !== w ? (r = [r, s, o, a, l, h, u, g], i = r) : (et = i, i = w)) : (et = i, i = w)) : (et = i, i = w)) : (et = i, i = w)) : (et = i, i = w)) : (et = i, i = w)) : (et = i, i = w)) : (et = i, i = w), i !== w && (tt = e, i = X(i)), e = i }

            function g() { var e, i, r, s; return e = et, i = et, 66 === t.charCodeAt(et) ? (r = G, et++) : (r = w, 0 === st && n(Y)), r !== w ? (114 === t.charCodeAt(et) ? (s = Z, et++) : (s = w, 0 === st && n(K)), s === w && (s = null), s !== w ? (r = [r, s], i = r) : (et = i, i = w)) : (et = i, i = w), i === w && (i = et, 67 === t.charCodeAt(et) ? (r = $, et++) : (r = w, 0 === st && n(J)), r !== w ? (108 === t.charCodeAt(et) ? (s = Q, et++) : (s = w, 0 === st && n(ee)), s === w && (s = null), s !== w ? (r = [r, s], i = r) : (et = i, i = w)) : (et = i, i = w), i === w && (te.test(t.charAt(et)) ? (i = t.charAt(et), et++) : (i = w, 0 === st && n(ie)))), i !== w && (tt = e, i = re(i)), e = i }

            function c() { var e, i; return e = et, ne.test(t.charAt(et)) ? (i = t.charAt(et), et++) : (i = w, 0 === st && n(se)), i !== w && (tt = e, i = M(i)), e = i }

            function v() { var e, i; return e = et, 42 === t.charCodeAt(et) ? (i = oe, et++) : (i = w, 0 === st && n(ae)), i !== w && (tt = e, i = le(i)), e = i }

            function d() { var e, i, r, s; return e = et, i = et, he.test(t.charAt(et)) ? (r = t.charAt(et), et++) : (r = w, 0 === st && n(ue)), r !== w ? (ge.test(t.charAt(et)) ? (s = t.charAt(et), et++) : (s = w, 0 === st && n(ce)), s === w && (s = null), s !== w ? (r = [r, s], i = r) : (et = i, i = w)) : (et = i, i = w), i !== w && (tt = e, i = ve(i)), e = i }

            function f() { var e, i, r, s, o; return e = et, i = et, 37 === t.charCodeAt(et) ? (r = de, et++) : (r = w, 0 === st && n(fe)), r !== w ? (pe.test(t.charAt(et)) ? (s = t.charAt(et), et++) : (s = w, 0 === st && n(ye)), s !== w ? (me.test(t.charAt(et)) ? (o = t.charAt(et), et++) : (o = w, 0 === st && n(be)), o !== w ? (r = [r, s, o], i = r) : (et = i, i = w)) : (et = i, i = w)) : (et = i, i = w), i === w && (me.test(t.charAt(et)) ? (i = t.charAt(et), et++) : (i = w, 0 === st && n(be))), i !== w && (tt = e, i = Ae(i)), e = i }

            function p() { var e, i, r, s, o, a, l; return e = et, i = et, 64 === t.charCodeAt(et) ? (r = xe, et++) : (r = w, 0 === st && n(ke)), r !== w ? (64 === t.charCodeAt(et) ? (s = xe, et++) : (s = w, 0 === st && n(ke)), s === w && (s = et, t.substr(et, 2) === Ce ? (o = Ce, et += 2) : (o = w, 0 === st && n(Re)), o !== w ? (we.test(t.charAt(et)) ? (a = t.charAt(et), et++) : (a = w, 0 === st && n(Se)), a !== w ? (o = [o, a], s = o) : (et = s, s = w)) : (et = s, s = w), s === w && (s = et, t.substr(et, 2) === Ve ? (o = Ve, et += 2) : (o = w, 0 === st && n(Ne)), o !== w ? (we.test(t.charAt(et)) ? (a = t.charAt(et), et++) : (a = w, 0 === st && n(Se)), a !== w ? (o = [o, a], s = o) : (et = s, s = w)) : (et = s, s = w), s === w && (s = et, t.substr(et, 2) === Le ? (o = Le, et += 2) : (o = w, 0 === st && n(Be)), o !== w ? (Te.test(t.charAt(et)) ? (a = t.charAt(et), et++) : (a = w, 0 === st && n(Ie)), a !== w ? (o = [o, a], s = o) : (et = s, s = w)) : (et = s, s = w), s === w && (s = et, t.substr(et, 2) === He ? (o = He, et += 2) : (o = w, 0 === st && n(Me)), o !== w ? (pe.test(t.charAt(et)) ? (a = t.charAt(et), et++) : (a = w, 0 === st && n(ye)), a !== w ? (me.test(t.charAt(et)) ? (l = t.charAt(et), et++) : (l = w, 0 === st && n(be)), l === w && (l = null), l !== w ? (o = [o, a, l], s = o) : (et = s, s = w)) : (et = s, s = w)) : (et = s, s = w), s === w && (s = et, t.substr(et, 2) === Pe ? (o = Pe, et += 2) : (o = w, 0 === st && n(Oe)), o !== w ? (pe.test(t.charAt(et)) ? (a = t.charAt(et), et++) : (a = w, 0 === st && n(ye)), a !== w ? (me.test(t.charAt(et)) ? (l = t.charAt(et), et++) : (l = w, 0 === st && n(be)), l === w && (l = null), l !== w ? (o = [o, a, l], s = o) : (et = s, s = w)) : (et = s, s = w)) : (et = s, s = w)))))), s === w && (s = null), s !== w ? (r = [r, s], i = r) : (et = i, i = w)) : (et = i, i = w), i !== w && (tt = e, i = We(i)), e = i }

            function y() { var e, t; return e = et, t = m(), t === w && (t = b()), t !== w && (tt = e, t = Fe(t)), e = t }

            function m() { var e, i, r, s, o, a; return e = et, i = et, 43 === t.charCodeAt(et) ? (r = Ee, et++) : (r = w, 0 === st && n(De)), r !== w ? (43 === t.charCodeAt(et) ? (s = Ee, et++) : (s = w, 0 === st && n(De)), s === w && (s = et, pe.test(t.charAt(et)) ? (o = t.charAt(et), et++) : (o = w, 0 === st && n(ye)), o !== w ? (me.test(t.charAt(et)) ? (a = t.charAt(et), et++) : (a = w, 0 === st && n(be)), a === w && (a = null), a !== w ? (o = [o, a], s = o) : (et = s, s = w)) : (et = s, s = w)), s === w && (s = null), s !== w ? (r = [r, s], i = r) : (et = i, i = w)) : (et = i, i = w), i !== w && (tt = e, i = ze(i)), e = i }

            function b() { var e, i, r, s, o, a; return e = et, i = et, 45 === t.charCodeAt(et) ? (r = _e, et++) : (r = w, 0 === st && n(qe)), r !== w ? (45 === t.charCodeAt(et) ? (s = _e, et++) : (s = w, 0 === st && n(qe)), s === w && (s = et, pe.test(t.charAt(et)) ? (o = t.charAt(et), et++) : (o = w, 0 === st && n(ye)), o !== w ? (me.test(t.charAt(et)) ? (a = t.charAt(et), et++) : (a = w, 0 === st && n(be)), a === w && (a = null), a !== w ? (o = [o, a], s = o) : (et = s, s = w)) : (et = s, s = w)), s === w && (s = null), s !== w ? (r = [r, s], i = r) : (et = i, i = w)) : (et = i, i = w), i !== w && (tt = e, i = Ue(i)), e = i }

            function A() { var e, i, r, s; return e = et, i = et, 72 === t.charCodeAt(et) ? (r = je, et++) : (r = w, 0 === st && n(Xe)), r !== w ? (me.test(t.charAt(et)) ? (s = t.charAt(et), et++) : (s = w, 0 === st && n(be)), s === w && (s = null), s !== w ? (r = [r, s], i = r) : (et = i, i = w)) : (et = i, i = w), i !== w && (tt = e, i = Ge(i)), e = i }

            function x() { var e, i, r, s, o, a, l; if (e = et, i = et, 58 === t.charCodeAt(et) ? (r = Ye, et++) : (r = w, 0 === st && n(Ze)), r !== w) { if (s = et, pe.test(t.charAt(et)) ? (o = t.charAt(et), et++) : (o = w, 0 === st && n(ye)), o !== w) { for (a = [], me.test(t.charAt(et)) ? (l = t.charAt(et), et++) : (l = w, 0 === st && n(be)); l !== w;) a.push(l), me.test(t.charAt(et)) ? (l = t.charAt(et), et++) : (l = w, 0 === st && n(be));
                        a !== w ? (o = [o, a], s = o) : (et = s, s = w) } else et = s, s = w;
                    s === w && (Ke.test(t.charAt(et)) ? (s = t.charAt(et), et++) : (s = w, 0 === st && n($e))), s !== w ? (r = [r, s], i = r) : (et = i, i = w) } else et = i, i = w; return i !== w && (tt = e, i = Je(i)), e = i }

            function k() { var e, i, r, s, o; return e = et, i = et, pe.test(t.charAt(et)) ? (r = t.charAt(et), et++) : (r = w, 0 === st && n(ye)), r !== w ? (me.test(t.charAt(et)) ? (s = t.charAt(et), et++) : (s = w, 0 === st && n(be)), s === w && (s = null), s !== w ? (me.test(t.charAt(et)) ? (o = t.charAt(et), et++) : (o = w, 0 === st && n(be)), o === w && (o = null), o !== w ? (r = [r, s, o], i = r) : (et = i, i = w)) : (et = i, i = w)) : (et = i, i = w), i !== w && (tt = e, i = Qe(i)), e = i }
            var C, R = arguments.length > 1 ? arguments[1] : {},
                w = {},
                S = { chain: o },
                V = o,
                N = function(e) { for (var t = [], i = [], r = 0; r < e[1].length; r++) t.push(e[1][r]); for (var r = 0; r < e[2].length; r++) { var n = e[2][r][0] ? e[2][r][0] : "-";
                        i.push({ bond: n, id: e[2][r][1] }) } for (var r = 0; r < e[3].length; r++) t.push(e[3][r]); for (var r = 0; r < e[6].length; r++) t.push(e[6][r]); return { atom: e[0], isBracket: !!e[0].element, branches: t, branchCount: t.length, ringbonds: i, ringbondCount: i.length, bond: e[4] ? e[4] : "-", next: e[5], hasNext: !!e[5] } },
                L = "(",
                B = { type: "literal", value: "(", description: '"("' },
                T = ")",
                I = { type: "literal", value: ")", description: '")"' },
                H = function(e) { var t = e[1] ? e[1] : "-"; return e[2].branchBond = t, e[2] },
                M = function(e) { return e },
                P = /^[\-=#$:\/\\.]/,
                O = { type: "class", value: "[-=#$:/\\\\.]", description: "[-=#$:/\\\\.]" },
                W = function(e) { return e },
                F = "[",
                E = { type: "literal", value: "[", description: '"["' },
                D = "se",
                z = { type: "literal", value: "se", description: '"se"' },
                _ = "as",
                q = { type: "literal", value: "as", description: '"as"' },
                U = "]",
                j = { type: "literal", value: "]", description: '"]"' },
                X = function(e) { return { isotope: e[1], element: e[2], chirality: e[3], hcount: e[4], charge: e[5], class: e[6] } },
                G = "B",
                Y = { type: "literal", value: "B", description: '"B"' },
                Z = "r",
                K = { type: "literal", value: "r", description: '"r"' },
                $ = "C",
                J = { type: "literal", value: "C", description: '"C"' },
                Q = "l",
                ee = { type: "literal", value: "l", description: '"l"' },
                te = /^[NOPSFI]/,
                ie = { type: "class", value: "[NOPSFI]", description: "[NOPSFI]" },
                re = function(e) { return e.length > 1 ? e.join("") : e },
                ne = /^[bcnops]/,
                se = {
                    type: "class",
                    value: "[bcnops]",
                    description: "[bcnops]"
                },
                oe = "*",
                ae = { type: "literal", value: "*", description: '"*"' },
                le = function(e) { return e },
                he = /^[A-Z]/,
                ue = { type: "class", value: "[A-Z]", description: "[A-Z]" },
                ge = /^[a-z]/,
                ce = { type: "class", value: "[a-z]", description: "[a-z]" },
                ve = function(e) { return e.join("") },
                de = "%",
                fe = { type: "literal", value: "%", description: '"%"' },
                pe = /^[1-9]/,
                ye = { type: "class", value: "[1-9]", description: "[1-9]" },
                me = /^[0-9]/,
                be = { type: "class", value: "[0-9]", description: "[0-9]" },
                Ae = function(e) { return 1 == e.length ? Number(e) : Number(e.join("").replace("%", "")) },
                xe = "@",
                ke = { type: "literal", value: "@", description: '"@"' },
                Ce = "TH",
                Re = { type: "literal", value: "TH", description: '"TH"' },
                we = /^[12]/,
                Se = { type: "class", value: "[12]", description: "[12]" },
                Ve = "AL",
                Ne = { type: "literal", value: "AL", description: '"AL"' },
                Le = "SP",
                Be = { type: "literal", value: "SP", description: '"SP"' },
                Te = /^[1-3]/,
                Ie = { type: "class", value: "[1-3]", description: "[1-3]" },
                He = "TB",
                Me = { type: "literal", value: "TB", description: '"TB"' },
                Pe = "OH",
                Oe = { type: "literal", value: "OH", description: '"OH"' },
                We = function(e) { return e[1] ? "@" == e[1] ? "@@" : e[1].join("").replace(",", "") : "@" },
                Fe = function(e) { return e },
                Ee = "+",
                De = { type: "literal", value: "+", description: '"+"' },
                ze = function(e) { return e[1] ? "+" != e[1] ? Number(e[1].join("")) : 2 : 1 },
                _e = "-",
                qe = { type: "literal", value: "-", description: '"-"' },
                Ue = function(e) { return e[1] ? "-" != e[1] ? -Number(e[1].join("")) : -2 : -1 },
                je = "H",
                Xe = { type: "literal", value: "H", description: '"H"' },
                Ge = function(e) { return e[1] ? Number(e[1]) : 1 },
                Ye = ":",
                Ze = { type: "literal", value: ":", description: '":"' },
                Ke = /^[0]/,
                $e = { type: "class", value: "[0]", description: "[0]" },
                Je = function(e) { return Number(e[1][0] + e[1][1].join("")) },
                Qe = function(e) { return Number(e.join("")) },
                et = 0,
                tt = 0,
                it = [{ line: 1, column: 1, seenCR: !1 }],
                rt = 0,
                nt = [],
                st = 0;
            if ("startRule" in R) { if (!(R.startRule in S)) throw new Error("Can't start parsing from rule \"" + R.startRule + '".');
                V = S[R.startRule] }
            if ((C = V()) !== w && et === t.length) return C;
            throw C !== w && et < t.length && n({ type: "end", description: "end of input" }), s(null, nt, rt < t.length ? t.charAt(rt) : null, rt < t.length ? r(rt, rt + 1) : r(rt, rt))
        }
        return function(e, t) {
            function i() { this.constructor = e }
            i.prototype = t.prototype, e.prototype = new i }(e, Error), { SyntaxError: e, parse: t }
    }(),
    SmilesDrawer = function() {
        function e(t) { _classCallCheck(this, e), this.ringIdCounter = 0, this.ringConnectionIdCounter = 0, this.canvasWrapper = null, this.direction = 1, this.totalOverlapScore = 0, this.maxBonds = { c: 4, C: 4, n: 3, N: 3, o: 2, O: 2 }, this.defaultOptions = { bondLength: 16, shortBondLength: 9, bondSpacing: 4, atomVisualization: "default", allowFlips: !1, isomeric: !1, debug: !1, themes: { dark: { C: "#fff", O: "#e74c3c", N: "#3498db", F: "#27ae60", CL: "#16a085", BR: "#d35400", I: "#8e44ad", P: "#d35400", S: "#f1c40f", B: "#e67e22", SI: "#e67e22", H: "#252525", BACKGROUND: "#141414" }, light: { C: "#222", O: "#e74c3c", N: "#3498db", F: "#27ae60", CL: "#16a085", BR: "#d35400", I: "#8e44ad", P: "#d35400", S: "#f1c40f", B: "#e67e22", SI: "#e67e22", H: "#d5d5d5", BACKGROUND: "#fff" } } }, this.opts = this.extend(!0, this.defaultOptions, t), this.theme = this.opts.themes.dark }
        return _createClass(e, [{ key: "extend", value: function() { var e = this,
                    t = {},
                    i = !1,
                    r = 0,
                    n = arguments.length; "[object Boolean]" === Object.prototype.toString.call(arguments[0]) && (i = arguments[0], r++); for (; r < n; r++) { var s = arguments[r];! function(r) { for (var n in r) Object.prototype.hasOwnProperty.call(r, n) && (i && "[object Object]" === Object.prototype.toString.call(r[n]) ? t[n] = e.extend(!0, t[n], r[n]) : t[n] = r[n]) }(s) } return t } }, { key: "draw", value: function(e, t) { var i = arguments.length > 2 && void 0 !== arguments[2] ? arguments[2] : "light",
                    r = arguments.length > 3 && void 0 !== arguments[3] && arguments[3]; if (this.data = e, this.canvasWrapper = new CanvasWrapper(t, this.opts.themes[i], this.opts.bondLength, this.opts.bondSpacing), this.ringIdCounter = 0, this.ringConnectionIdCounter = 0, this.vertices = [], this.edges = [], this.rings = [], this.ringConnections = [], this.originalRings = [], this.originalRingConnections = [], this.bridgedRing = !1, this.initGraph(e), this.initRings(), this.opts.isomeric && this.annotateChirality(), !r) { this.position(), this.restoreRingInformation(); var n = this.getOverlapScore();
                    this.totalOverlapScore = this.getOverlapScore().total; for (var s = 0; s < this.edges.length; s++) { var o = this.edges[s]; if (this.isEdgeRotatable(o)) { var a = this.getTreeDepth(o.sourceId, o.targetId),
                                l = this.getTreeDepth(o.targetId, o.sourceId),
                                h = o.targetId,
                                u = o.sourceId;
                            a > l && (h = o.sourceId, u = o.targetId, l); if (this.getSubtreeOverlapScore(u, h, n.vertexScores).value > 1) { var g = this.vertices[h],
                                    c = this.vertices[u],
                                    v = c.getNeighbours(h); if (1 === v.length) { var d = this.vertices[v[0]],
                                        f = d.position.getRotateAwayFromAngle(g.position, c.position, MathHelper.toRad(120));
                                    this.rotateSubtree(d.id, c.id, f, c.position); var p = this.getOverlapScore().total;
                                    p > this.totalOverlapScore ? this.rotateSubtree(d.id, c.id, -f, c.position) : this.totalOverlapScore = p } else if (2 == v.length) { if (c.value.rings.length + g.value.rings.length > 0) continue; var y = this.vertices[v[0]],
                                        m = this.vertices[v[1]],
                                        b = y.position.getRotateAwayFromAngle(g.position, c.position, MathHelper.toRad(120)),
                                        A = m.position.getRotateAwayFromAngle(g.position, c.position, MathHelper.toRad(120));
                                    this.rotateSubtree(y.id, c.id, b, c.position), this.rotateSubtree(m.id, c.id, A, c.position); var x = this.getOverlapScore().total;
                                    x > this.totalOverlapScore ? (this.rotateSubtree(y.id, c.id, -b, c.position), this.rotateSubtree(m.id, c.id, -A, c.position)) : this.totalOverlapScore = x }
                                n = this.getOverlapScore() } } }
                    this.resolveSecondaryOverlaps(n.scores), this.canvasWrapper.scale(this.vertices), this.drawEdges(this.opts.debug), this.drawVertices(this.opts.debug), this.canvasWrapper.reset() } } }, { key: "edgeRingCount", value: function(e) { var t = this.edges[e],
                    i = this.vertices[t.sourceId],
                    r = this.vertices[t.targetId]; return Math.min(i.value.rings.length, r.value.rings.length) } }, { key: "getBridgedRings", value: function() { for (var e = [], t = 0; t < this.rings.length; t++) this.rings[t].isBridged && e.push(this.rings[t]); return e } }, { key: "getFusedRings", value: function() { for (var e = [], t = 0; t < this.rings.length; t++) this.rings[t].isFused && e.push(this.rings[t]); return e } }, { key: "getSpiros", value: function() { for (var e = [], t = 0; t < this.rings.length; t++) this.rings[t].isSpiro && e.push(this.rings[t]); return e } }, { key: "printRingInfo", value: function() { for (var e = "", t = 0; t < this.rings.length; t++) { var i = this.rings[t];
                    e += i.id + ";", e += i.members.length + ";", e += i.neighbours.length + ";", e += i.isSpiro ? "true;" : "false;", e += i.isFused ? "true;" : "false;", e += i.isBridged ? "true;" : "false;", e += i.rings.length + ";", e += i.insiders.length, e += "\n" } return e } }, { key: "getTotalOverlapScore", value: function() { return this.totalOverlapScore } }, { key: "getRingCount", value: function() { return this.rings.length } }, { key: "hasBridgedRing", value: function() { return this.bridgedRing } }, { key: "getHeavyAtomCount", value: function() { for (var e = 0, t = 0; t < this.vertices.length; t++) "h" !== this.vertices[t].value.element.toLowerCase() && e++; return e } }, { key: "initGraph", value: function(e) { var t = arguments.length > 1 && void 0 !== arguments[1] ? arguments[1] : 0,
                    i = arguments.length > 2 && void 0 !== arguments[2] ? arguments[2] : null,
                    r = arguments.length > 3 && void 0 !== arguments[3] && arguments[3],
                    n = new Atom(e.atom.element ? e.atom.element : e.atom, e.bond);
                n.branchBond = e.branchBond, n.ringbonds = e.ringbonds, n.bracket = e.atom.element ? e.atom : null, n.setOrder(i, t); var s = new Vertex(n),
                    o = this.vertices[i]; if (s.parentVertexId = i, this.addVertex(s), null != i) { this.vertices[i].children.push(s.id), this.vertices[i].spanningTreeChildren.push(s.id); var a = new Edge(i, s.id, 1);
                    a.bondType = r ? s.value.branchBond : this.vertices[i].value.bondType; var l = this.addEdge(a);
                    s.edges.push(l), o.edges.push(l) } if (n.bracket && this.opts.isomeric)
                    for (var h = 0; h < n.bracket.hcount; h++) this.opts.isomeric && this.initGraph({ atom: { element: "H", bond: "-" }, ringbonds: [] }, h + 1, s.id); var u = e.ringbondCount + 1;
                n.bracket && (u += n.bracket.hcount); for (var g = 0; g < e.branchCount; g++) this.initGraph(e.branches[g], g + u, s.id, !0);
                e.hasNext && this.initGraph(e.next, e.branchCount + u, s.id) } }, { key: "getRingbondType", value: function(e, t) { if (e.value.getRingbondCount() < 1 || t.value.getRingbondCount() < 1) return null; for (var i = 0; i < e.value.ringbonds.length; i++)
                    for (var r = 0; r < t.value.ringbonds.length; r++)
                        if (e.value.ringbonds[i].id === t.value.ringbonds[r].id) return "-" === e.value.ringbonds[i].bondType ? t.value.ringbonds[r].bond : e.value.ringbonds[i].bond;
                return null } }, { key: "initRings", value: function() { for (var e = this, t = {}, i = this.vertices.length - 1; i >= 0; i--) { var r = this.vertices[i]; if (0 !== r.value.ringbonds.length)
                        for (var n = 0; n < r.value.ringbonds.length; n++) { var s = r.value.ringbonds[n].id; if (void 0 === t[s]) t[s] = r.id;
                            else { var o = t[s],
                                    a = r.id,
                                    l = e.addEdge(new Edge(a, o, 1)),
                                    h = e.vertices[a],
                                    u = e.vertices[o];
                                h.children.push(o), u.children.push(a), h.edges.push(l), u.edges.push(l); var g = new Ring(s, a, o);
                                e.addRing(g); for (var c = e.getRingVertices(g.sourceId, g.targetId), v = 0; v < c.length; v++) g.members.push(c[v]), e.vertices[c[v]].value.rings.push(g.id);
                                t[s] = void 0; var d = u.value.bracket ? u.value.bracket.hcount : 0,
                                    f = h.value.bracket ? h.value.bracket.hcount : 0;
                                u.value.setOrder(a, n + 1 + f), h.value.setOrder(o, n + 1 + d) } } } for (var p = 0; p < this.rings.length - 1; p++)
                    for (var y = p + 1; y < this.rings.length; y++) { var m = this.rings[p],
                            b = this.rings[y],
                            A = new RingConnection(m, b);
                        A.vertices.length > 0 && this.addRingConnection(A) }
                for (var x = 0; x < this.rings.length; x++) { var k = this.rings[x];
                    k.neighbours = RingConnection.getNeighbours(this.ringConnections, k.id) } for (this.backupRingInformation(); this.rings.length > 0;) { for (var C = -1, R = 0; R < this.rings.length; R++) { var w = this.rings[R];
                        this.isPartOfBridgedRing(w.id) && (C = w.id) } if (-1 === C) break; var S = this.getRing(C),
                        V = this.getBridgedRingRings(S.id);
                    this.bridgedRing = !0, this.createBridgedRing(V, S.sourceId); for (var N = 0; N < V.length; N++) this.removeRing(V[N]) } } }, { key: "getBridgedRingRings", value: function(e) { var t = new Array,
                    i = this; return function e(r) { var n = i.getRing(r);
                    t.push(r); for (var s = 0; s < n.neighbours.length; s++) { var o = n.neighbours[s]; - 1 === t.indexOf(o) && o !== r && RingConnection.isBridge(i.ringConnections, i.vertices, r, o) && e(o) } }(e), ArrayHelper.unique(t) } }, { key: "isPartOfBridgedRing", value: function(e) { for (var t = 0; t < this.ringConnections.length; t++)
                    if (this.ringConnections[t].rings.contains(e) && this.ringConnections[t].isBridge(this.vertices)) return !0;
                return !1 } }, { key: "createBridgedRing", value: function(e, t) { for (var i = new Array, r = new Array, n = new Array, s = (new Array, 0); s < e.length; s++) { for (var o = this.getRing(e[s]), a = 0; a < o.members.length; a++) r.push(o.members[a]); for (var l = 0; l < o.neighbours.length; l++) n.push(o.neighbours[l]) }
                r = ArrayHelper.unique(r); for (var h = new Array, u = 0; u < r.length; u++) { var g = this.vertices[r[u]],
                        c = ArrayHelper.intersection(e, g.value.rings);
                    1 == g.value.rings.length || 1 == c.length ? i.push(g.id) : h.push(g.id) } for (var v = new Array, d = new Array, f = 0; f < h.length; f++) { for (var p = this.vertices[h[f]], y = !1, m = 0; m < p.edges.length; m++) 1 == this.edgeRingCount(p.edges[m]) && (y = !0);
                    y ? (p.value.isBridgeNode = !0, v.push(p.id)) : (p.value.isBridge = !0, d.push(p.id)) } var b = ArrayHelper.merge(i, v);
                n = ArrayHelper.unique(n), n = ArrayHelper.removeAll(n, e); for (var A = this.vertices[t], x = A.getNeighbours(), k = null, C = 0; C < x.length; C++) { var R = x[C]; - 1 !== b.indexOf(R) && (k = R) } var w = new Ring(-1, t, k);
                w.isBridged = !0, w.members = b, w.neighbours = n, w.insiders = d; for (var S = 0; S < e.length; S++) w.rings.push(this.getRing(e[S]).clone());
                this.addRing(w), this.vertices[t].value.anchoredRings.push(w.id); for (var V = 0; V < d.length; V++) { var N = this.vertices[d[V]];
                    N.value.rings = new Array, N.value.anchoredRings = new Array, N.value.bridgedRing = w.id } for (var L = 0; L < b.length; L++) { var B = this.vertices[b[L]];
                    B.value.rings = ArrayHelper.removeAll(B.value.rings, e), B.value.rings.push(w.id) } for (var T = 0; T < e.length; T++)
                    for (var I = T + 1; I < e.length; I++) this.removeRingConnectionsBetween(e[T], e[I]); for (var H = 0; H < n.length; H++) { for (var M = this.getRingConnections(n[H], e), P = 0; P < M.length; P++) this.getRingConnection(M[P]).updateOther(w.id, n[H]);
                    this.getRing(n[H]).neighbours.push(w.id) } return w } }, { key: "getRingVertices", value: function(e, t) { for (var i = this.dijkstra(e, t), r = [], n = [], s = t; null != s;) r.push(s), s = i[s]; for (var o = r.length - 1; o >= 0; o--) n.push(r[o]); return n } }, { key: "dijkstra", value: function(e, t) { for (var i = new Array(this.vertices.length), r = new Array(this.vertices.length), n = new Array(this.vertices.length), s = new Array(this.vertices.length), o = 0; o < this.vertices.length; o++) r[o] = o === e ? 0 : Number.MAX_VALUE, i[o] = null, n[o] = !1, s[o] = this.vertices[o].getNeighbours(); for (; ArrayHelper.count(n, !1) > 0;) { var a = this.getMinDist(r, n); if (a == t) return i;
                    n[a] = !0; for (var l = 0; l < s[a].length; l++) { var h = s[a][l],
                            u = r[a] + this.getEdgeWeight(a, h);
                        a == e && h == t || a == t && h == e || u < r[h] && (r[h] = u, i[h] = a) } } } }, { key: "getMinDist", value: function(e, t) { for (var i = Number.MAX_VALUE, r = null, n = 0; n < e.length; n++) t[n] || e[n] < i && (r = n, i = e[r]); return r } }, { key: "areVerticesInSameRing", value: function(e, t) { for (var i = 0; i < e.value.rings.length; i++)
                    for (var r = 0; r < t.value.rings.length; r++)
                        if (e.value.rings[i] == t.value.rings[r]) return !0;
                return !1 } }, { key: "getCommonRings", value: function(e, t) { for (var i = [], r = 0; r < e.value.rings.length; r++)
                    for (var n = 0; n < t.value.rings.length; n++) e.value.rings[r] == t.value.rings[n] && i.push(e.value.rings[r]); return i } }, { key: "getSmallestCommonRing", value: function(e, t) { for (var i = this.getCommonRings(e, t), r = Number.MAX_VALUE, n = null, s = 0; s < i.length; s++) { var o = this.getRing(i[s]).getSize();
                    o < r && (r = o, n = this.getRing(i[s])) } return n } }, { key: "getLargestCommonRing", value: function(e, t) { for (var i = this.getCommonRings(e, t), r = 0, n = null, s = 0; s < i.length; s++) { var o = this.getRing(i[s]).getSize();
                    o > r && (r = o, n = this.getRing(i[s])) } return n } }, { key: "getLargestOrAromaticCommonRing", value: function(e, t) { for (var i = this.getCommonRings(e, t), r = 0, n = null, s = 0; s < i.length; s++) { var o = this.getRing(i[s]),
                        a = o.getSize(); if (o.isBenzeneLike(this.vertices)) return o;
                    a > r && (r = a, n = o) } return n } }, { key: "getVerticesAt", value: function(e, t, i) { for (var r = new Array, n = 0; n < this.vertices.length; n++) { var s = this.vertices[n]; if (s.id !== i && s.positioned) { e.distance(s.position) <= t && r.push(s.id) } } return r } }, { key: "getClosestVertex", value: function(e) { for (var t = 99999, i = null, r = 0; r < this.vertices.length; r++) { var n = this.vertices[r]; if (n.id !== e.id) { var s = e.position.distanceSq(n.position);
                        s < t && (t = s, i = n) } } return i } }, { key: "getClosestEndpointVertex", value: function(e) { for (var t = 99999, i = null, r = 0; r < this.vertices.length; r++) { var n = this.vertices[r]; if (!(n.id === e.id || n.getNeighbours().length > 1)) { var s = e.position.distanceSq(n.position);
                        s < t && (t = s, i = n) } } return i } }, { key: "getBranch", value: function(e, t) { var i = new Array,
                    r = new Array,
                    n = this; return i.push(e),
                    function e(t, s) { for (var o = n.vertices[t], a = 0; a < o.value.rings.length; a++) r.push(o.value.rings[a]); for (var l = 0; l < o.children.length; l++) { var h = o.children[l];
                            h === s || ArrayHelper.contains(i, { value: h }) || (i.push(h), e(h, t)) } var u = o.parentVertexId;
                        u === s || null === u || ArrayHelper.contains(i, { value: u }) || (i.push(u), e(u, t)) }(e, t), { vertices: i, rings: ArrayHelper.unique(r) } } }, { key: "addVertex", value: function(e) { return e.id = this.vertices.length, this.vertices.push(e), e.id } }, { key: "addEdge", value: function(e) { return e.id = this.edges.length, this.edges.push(e), e.id } }, { key: "addRing", value: function(e) { return e.id = this.ringIdCounter++, this.rings.push(e), e.id } }, { key: "removeRing", value: function(e) { this.rings = this.rings.filter(function(t) { return t.id !== e }), this.ringConnections = this.ringConnections.filter(function(t) { return t.rings.first !== e && t.rings.second !== e }); for (var t = 0; t < this.rings.length; t++) { var i = this.rings[t];
                    i.neighbours = i.neighbours.filter(function(t) { return t !== e }) } } }, { key: "getRing", value: function(e) { for (var t = 0; t < this.rings.length; t++)
                    if (this.rings[t].id == e) return this.rings[t] } }, { key: "addRingConnection", value: function(e) { return e.id = this.ringConnectionIdCounter++, this.ringConnections.push(e), e.id } }, { key: "removeRingConnection", value: function(e) { this.ringConnections = this.ringConnections.filter(function(t) { return t.id !== e }) } }, { key: "removeRingConnectionsBetween", value: function(e, t) { for (var i = new Array, r = 0; r < this.ringConnections.length; r++) { var n = this.ringConnections[r];
                    (n.rings.first === e && n.rings.second === t || n.rings.first === t && n.rings.second === e) && i.push(n.id) } for (var s = 0; s < i.length; s++) this.removeRingConnection(i[s]) } }, { key: "getRingConnection", value: function(e) { for (var t = 0; t < this.ringConnections.length; t++)
                    if (this.ringConnections[t].id == e) return this.ringConnections[t] } }, { key: "getRingConnections", value: function(e) { var t = arguments.length > 1 && void 0 !== arguments[1] ? arguments[1] : null,
                    i = new Array; if (null === t)
                    for (var r = 0; r < this.ringConnections.length; r++) { var n = this.ringConnections[r];
                        n.rings.first !== e && n.rings.second !== e || i.push(n.id) } else if (t.constructor !== Array)
                        for (var s = 0; s < this.ringConnections.length; s++) { var o = this.ringConnections[s];
                            (o.rings.first === e && o.rings.second === t || o.rings.first === t && o.rings.second === e) && i.push(o.id) } else
                            for (var a = 0; a < this.ringConnections.length; a++)
                                for (var l = 0; l < t.length; l++) { var h = t[l],
                                        u = this.ringConnections[a];
                                    (u.rings.first === e && u.rings.second === h || u.rings.first === h && u.rings.second === e) && i.push(u.id) }
                    return i } }, { key: "isRingConnection", value: function(e, t) { for (var i = 0; i < this.ringConnections.length; i++) { var r = this.ringConnections[i]; if (2 === r.vertices.length && (r.vertices[0] === e && r.vertices[1] === t || r.vertices[0] === t && r.vertices[1] === e)) return !0 } return !1 } }, { key: "getOverlapScore", value: function() { for (var e = 0, t = new Float32Array(this.vertices.length), i = 0; i < this.vertices.length; i++) t[i] = 0; for (var r = 0; r < this.vertices.length; r++)
                    for (var n = r + 1; n < this.vertices.length; n++) { var s = this.vertices[r],
                            o = this.vertices[n],
                            a = Vector2.subtract(s.position, o.position).length(); if (a < this.opts.bondLength) { var l = (this.opts.bondLength - a) / this.opts.bondLength;
                            e += l, t[r] += l, t[n] += l } }
                for (var h = [], u = 0; u < this.vertices.length; u++) h.push({ id: u, score: t[u] }); return h.sort(function(e, t) { return t.score - e.score }), { total: e, scores: h, vertexScores: t } } }, { key: "chooseSide", value: function(e, t, i) { for (var r = e.getNeighbours(t.id), n = t.getNeighbours(e.id), s = r.length, o = n.length, a = ArrayHelper.merge(r, n), l = [0, 0], h = 0; h < a.length; h++) { this.vertices[a[h]].position.sameSideAs(e.position, t.position, i[0]) ? l[0]++ : l[1]++ } for (var u = [0, 0], g = 0; g < this.vertices.length; g++) { this.vertices[g].position.sameSideAs(e.position, t.position, i[0]) ? u[0]++ : u[1]++ } return { totalSideCount: u, totalPosition: u[0] > u[1] ? 0 : 1, sideCount: l, position: l[0] > l[1] ? 0 : 1, anCount: s, bnCount: o } } }, { key: "areConnected", value: function(e, t) { for (var i = 0; i < this.edges.length; i++) { var r = this.edges[i]; if (r.sourceId === e && r.targetId === t || r.sourceId === t && r.targetId === e) return !0 } return !1 } }, { key: "getEdgeWeight", value: function(e, t) { for (var i = 0; i < this.edges.length; i++) { var r = this.edges[i]; if (r.sourceId == e && r.targetId == t || r.targetId == e && r.sourceId == t) return r.weight } return null } }, { key: "getEdge", value: function(e, t) { for (var i = 0; i < this.edges.length; i++) { var r = this.edges[i]; if (r.sourceId == e && r.targetId == t || r.targetId == e && r.sourceId == t) return r } return null } }, { key: "forceLayout", value: function(e, t, i, r) { for (var n = this.opts.bondLength, s = this.vertices[i], o = s.getNeighbours(), a = 0; a < o.length; a++) this.vertices[o[a]].positioned && e.push(o[a]); for (var l = e.length + r.rings.length, h = new Array(e.length), u = {}, g = new Array(l), c = new Array, v = 0; v < l; v++) { g[v] = new Array(l); for (var d = 0; d < l; d++) g[v][d] = 0 } for (var f = 0; f < e.length; f++) h[f] = this.vertices[e[f]].id, u[h[f]] = f; for (var p = 0; p < e.length - 1; p++)
                    for (var y = p; y < e.length; y++) { var m = this.getEdge(h[p], this.vertices[e[y]].id);
                        null !== m && (g[p][y] = n, g[y][p] = n, c.push([p, y])) }
                for (var b = 0; b < r.rings.length; b++)
                    for (var A = r.rings[b], x = e.length + b, k = 0; k < A.members.length; k++) { var C = u[A.members[k]],
                            R = MathHelper.polyCircumradius(n, A.getSize());
                        g[C][x] = R, g[x][C] = R }
                for (var w = 0; w < c.length; w++) { for (var S = 0; S < l; S++) g[S].push(0);
                    g.push(new Array); for (var V = 0; V < l + c.length; V++) g[l + w].push(0) } for (var N = 0; N < r.rings.length; N++) { for (var L = r.rings[N], B = e.length + N, T = L.getSize(), I = 0; I < c.length; I++) { var H = c[I][0]; if (0 !== g[B][H]) { var M = MathHelper.apothem(g[B][H], T);
                            g[B][l + I] = M, g[l + I][B] = M } } for (var P = 0; P < r.rings.length; P++) { var O = r.rings[P]; if (O.id !== L.id) { if (0 !== ArrayHelper.intersection(L.members, O.members).length) { var W = e.length + P,
                                    F = O.getSize(),
                                    E = MathHelper.apothemFromSideLength(n, T) + MathHelper.apothemFromSideLength(n, F);
                                g[B][W] = E } } } }
                l += c.length; for (var D = l - c.length, z = new Array(l), _ = new Array(l), q = new Array(l), U = new Array(l), j = new Array(l), X = new Array(l), G = 0; G < l; G++) U[G] = G >= e.length && G < D, X[G] = G < e.length ? this.vertices[h[G]].value.originalRings.length : 1, U[G] ? j[G] = r.rings[G - e.length].members.length : j[G] = 1; for (var Y = 0; Y < l; Y++)
                    if (z[Y] = new Vector2, _[Y] = new Vector2(t.x + Math.random() * n, t.y + Math.random() * n), q[Y] = !1, !(Y >= e.length)) { var Z = this.vertices[h[Y]];
                        _[Y] = Z.position.clone(), 0 === Z.position.x && Z.position.y, Z.positioned && 2 === r.rings.length && (q[Y] = !0) }
                for (var K = n / 1.4, $ = n / 2, J = 2 * n, Q = 0; Q < 600; Q++) { for (var ee = 0; ee < l; ee++) z[ee].set(0, 0); for (var te = 0; te < c.length; te++) { var ie = D + te,
                            re = _[c[te][0]],
                            ne = _[c[te][1]];
                        _[ie] = Vector2.midpoint(re, ne) } for (var se = 0; se < l - 1; se++)
                        for (var oe = se + 1; oe < l; oe++)
                            if ((!(Q <= 250) || U[se] && U[oe]) && !(Q > 250 && U[se] && U[oe] || r.rings.length < 3 && (U[se] || U[oe]))) { var ae = _[oe].x - _[se].x,
                                    le = _[oe].y - _[se].y; if (0 !== ae && 0 !== le) { var he = ae * ae + le * le;
                                    he < .01 && (ae = .1 * Math.random() + .1, le = .1 * Math.random() + .1, he = ae * ae + le * le); var ue = Math.sqrt(he); if (!(ue > g[se][oe] && Q > 200)) { var ge = K * K / ue;
                                        Q <= 200 && (ge *= j[se] * j[oe]), Q > 250 && (U[se] || U[oe]) && (ge *= j[se] * j[oe]); var ce = ge * ae / ue,
                                            ve = ge * le / ue;
                                        q[se] || (z[se].x -= ce, z[se].y -= ve), q[oe] || (z[oe].x += ce, z[oe].y += ve) } } }
                    for (var de = 0; de < l - 1; de++)
                        for (var fe = de + 1; fe < l; fe++)
                            if (!(g[de][fe] <= 0) && (!(Q <= 250) || U[de] && U[fe]) && !(Q > 250 && U[de] && U[fe])) { var pe = _[fe].x - _[de].x,
                                    ye = _[fe].y - _[de].y; if (0 !== pe && 0 !== ye) { var me = pe * pe + ye * ye;
                                    me < .01 && (pe = .1 * Math.random() + .1, ye = .1 * Math.random() + .1, me = pe * pe + ye * ye); var be = Math.sqrt(me);
                                    be > J && (be = J, me = be * be); var Ae = (me - K * K) / K,
                                        xe = g[de][fe];
                                    Ae *= be / xe; var ke = Ae * pe / be,
                                        Ce = Ae * ye / be;
                                    q[de] || (z[de].x += ke, z[de].y += Ce), q[fe] || (z[fe].x -= ke, z[fe].y -= Ce) } }
                    for (var Re = 0; Re < c.length; Re++) { var we = D + Re,
                            Se = z[we],
                            Ve = c[Re][0],
                            Ne = c[Re][1];
                        z[Ve].x += Se.x, z[Ve].y += Se.y, z[Ne].x += Se.x, z[Ne].y += Se.y } for (var Le = 0; Le < l; Le++)
                        if (!q[Le]) { var Be = .005 * z[Le].x,
                                Te = .005 * z[Le].y;
                            Be > $ && (Be = $), Be < -$ && (Be = -$), Te > $ && (Te = $), Te < -$ && (Te = -$);
                            _[Le].x += Be, _[Le].y += Te }
                    if (Q > 200 && r.rings.length > 2)
                        for (var Ie = 0; Ie < r.rings.length; Ie++) { for (var He = r.rings[Ie], Me = new Vector2, Pe = 0; Pe < He.members.length; Pe++) { var Oe = _[u[He.members[Pe]]];
                                Me.x += Oe.x, Me.y += Oe.y }
                            Me.x /= He.members.length, Me.y /= He.members.length, _[e.length + Ie] = Me } } for (var We = 0; We < l; We++)
                    if (We < e.length) q[We] || (this.vertices[h[We]].position = _[We], this.vertices[h[We]].positioned = !0);
                    else if (We < e.length + r.rings.length) { var Fe = We - e.length;
                    r.rings[Fe].center = _[We] } for (var Ee = 0; Ee < e.length; Ee++)
                    for (var De = this.vertices[e[Ee]], ze = (this.vertices[De.parentVertexId], De.getNeighbours()), _e = 0; _e < ze.length; _e++) { var qe = this.vertices[ze[_e]];
                        qe.positioned || (t = this.getSubringCenter(r, De), this.createNextBond(qe, De, t)) }
                this.createRing(r, null, null, null, !0) } }, { key: "getSubringCenter", value: function(e, t) { for (var i = Number.MAX_VALUE, r = e.center, n = 0; n < e.rings.length; n++)
                    for (var s = e.rings[n], o = 0; o < s.members.length; o++) s.members[o] === t.id && i > s.members.length && (r = s.center, i = s.members.length); return r } }, { key: "drawEdges", value: function(e) { for (var t = this, i = this, r = 0; r < this.edges.length; r++) ! function(r) { var n = t.edges[r],
                        s = t.vertices[n.sourceId],
                        o = t.vertices[n.targetId],
                        a = s.value.element,
                        l = o.value.element,
                        h = s.position,
                        u = o.position,
                        g = t.getEdgeNormals(n),
                        c = ArrayHelper.clone(g); if (ArrayHelper.each(c, function(e) { e.multiply(10), e.add(h) }), "=" === n.bondType || "=" === t.getRingbondType(s, o)) { var v = t.areVerticesInSameRing(s, o),
                            d = t.chooseSide(s, o, c); if (v) { var f = t.getLargestOrAromaticCommonRing(s, o),
                                p = f.center;
                            ArrayHelper.each(g, function(e) { e.multiply(i.opts.bondSpacing) }); var y = null;
                            y = p.sameSideAs(s.position, o.position, Vector2.add(h, g[0])) ? new Line(Vector2.add(h, g[0]), Vector2.add(u, g[0]), a, l) : new Line(Vector2.add(h, g[1]), Vector2.add(u, g[1]), a, l), y.shorten(t.opts.bondLength - t.opts.shortBondLength), t.canvasWrapper.drawLine(y), t.canvasWrapper.drawLine(new Line(h, u, a, l)) } else if (n.center) { ArrayHelper.each(g, function(e) { e.multiply(i.opts.bondSpacing / 2) }); var m = new Line(Vector2.add(h, g[0]), Vector2.add(u, g[0]), a, l),
                                b = new Line(Vector2.add(h, g[1]), Vector2.add(u, g[1]), a, l);
                            m.shorten(t.opts.bondLength - t.opts.shortBondLength), b.shorten(t.opts.bondLength - t.opts.shortBondLength), t.canvasWrapper.drawLine(m), t.canvasWrapper.drawLine(b) } else if (0 == d.anCount && d.bnCount > 1 || 0 == d.bnCount && d.anCount > 1) { ArrayHelper.each(g, function(e) { e.multiply(i.opts.bondSpacing / 2) }); var A = new Line(Vector2.add(h, g[0]), Vector2.add(u, g[0]), a, l),
                                x = new Line(Vector2.add(h, g[1]), Vector2.add(u, g[1]), a, l);
                            t.canvasWrapper.drawLine(A), t.canvasWrapper.drawLine(x) } else if (d.sideCount[0] > d.sideCount[1]) { ArrayHelper.each(g, function(e) { e.multiply(i.opts.bondSpacing) }); var k = new Line(Vector2.add(h, g[0]), Vector2.add(u, g[0]), a, l);
                            k.shorten(t.opts.bondLength - t.opts.shortBondLength), t.canvasWrapper.drawLine(k), t.canvasWrapper.drawLine(new Line(h, u, a, l)) } else if (d.sideCount[0] < d.sideCount[1]) { ArrayHelper.each(g, function(e) { e.multiply(i.opts.bondSpacing) }); var C = new Line(Vector2.add(h, g[1]), Vector2.add(u, g[1]), a, l);
                            C.shorten(t.opts.bondLength - t.opts.shortBondLength), t.canvasWrapper.drawLine(C), t.canvasWrapper.drawLine(new Line(h, u, a, l)) } else if (d.totalSideCount[0] > d.totalSideCount[1]) { ArrayHelper.each(g, function(e) { e.multiply(i.opts.bondSpacing) }); var R = new Line(Vector2.add(h, g[0]), Vector2.add(u, g[0]), a, l);
                            R.shorten(t.opts.bondLength - t.opts.shortBondLength), t.canvasWrapper.drawLine(R), t.canvasWrapper.drawLine(new Line(h, u, a, l)) } else if (d.totalSideCount[0] <= d.totalSideCount[1]) { ArrayHelper.each(g, function(e) { e.multiply(i.opts.bondSpacing) }); var w = new Line(Vector2.add(h, g[1]), Vector2.add(u, g[1]), a, l);
                            w.shorten(t.opts.bondLength - t.opts.shortBondLength), t.canvasWrapper.drawLine(w), t.canvasWrapper.drawLine(new Line(h, u, a, l)) } } else if ("#" === n.bondType) { ArrayHelper.each(g, function(e) { e.multiply(i.opts.bondSpacing / 1.5) }); var S = new Line(Vector2.add(h, g[0]), Vector2.add(u, g[0]), a, l),
                            V = new Line(Vector2.add(h, g[1]), Vector2.add(u, g[1]), a, l);
                        S.shorten(t.opts.bondLength - t.opts.shortBondLength), V.shorten(t.opts.bondLength - t.opts.shortBondLength), t.canvasWrapper.drawLine(S), t.canvasWrapper.drawLine(V), t.canvasWrapper.drawLine(new Line(h, u, a, l)) } else { var N = s.value.bracket && s.value.bracket.chirality,
                            L = o.value.bracket && o.value.bracket.chirality; "up" === n.chiral ? t.canvasWrapper.drawWedge(new Line(h, u, a, l, N, L)) : "down" === n.chiral ? t.canvasWrapper.drawDashedWedge(new Line(h, u, a, l, N, L)) : t.canvasWrapper.drawLine(new Line(h, u, a, l, N, L)) } if (e) { var B = Vector2.midpoint(h, u);
                        t.canvasWrapper.drawDebugText(B.x, B.y, "e: " + r) } }(r); for (var r = 0; r < this.rings.length; r++) { var n = this.rings[r];
                    n.isAromatic(this.vertices) && this.canvasWrapper.drawAromaticityRing(n) } } }, { key: "drawVertices", value: function(e) { for (var t = 0; t < this.vertices.length; t++) { var i = this.vertices[t],
                        r = i.value,
                        n = 0,
                        s = this.getBondCount(i),
                        o = 1 == r.element.length ? r.element.toUpperCase() : r.element,
                        a = this.maxBonds[o] - s,
                        l = i.getTextDirection(this.vertices),
                        h = i.isTerminal(),
                        u = "c" === r.element.toLowerCase(); if (r.bracket && (a = r.bracket.hcount, n = r.bracket.charge), (!u || r.explicit || h) && ("default" === this.opts.atomVisualization ? this.canvasWrapper.drawText(i.position.x, i.position.y, o, a, l, h, n) : "balls" === this.opts.atomVisualization && this.canvasWrapper.drawBall(i.position.x, i.position.y, o)), e) { var g = "v: " + i.id + " " + ArrayHelper.print(r.ringbonds);
                        this.canvasWrapper.drawDebugText(i.position.x, i.position.y, g) } } if (this.opts.debug)
                    for (var c = 0; c < this.rings.length; c++) { var v = this.rings[c].center;
                        this.canvasWrapper.drawDebugPoint(v.x, v.y, "r: " + this.rings[c].id) } } }, { key: "position", value: function() { for (var e = this.vertices[0], t = 0; t < this.rings.length; t++)
                    if (this.rings[t].isBridged)
                        for (var i = 0; i < this.rings[t].members.length && (e = this.vertices[this.rings[t].members[i]], 1 !== e.value.originalRings.length); i++);
                this.createNextBond(e), this.resolvePrimaryOverlaps() } }, { key: "clearPositions", value: function() { this.vertexPositionsBackup = [], this.ringPositionsBackup = []; for (var e = 0; e < this.vertices.length; e++) { var t = this.vertices[e];
                    this.vertexPositionsBackup.push(t.position.clone()), t.positioned = !1, t.position = new Vector2 } for (var i = 0; i < this.rings.length; i++) { var r = this.rings[i];
                    this.ringPositionsBackup.push(r.center.clone()), r.positioned = !1, r.center = new Vector2 } } }, { key: "restorePositions", value: function() { for (var e = 0; e < this.vertexPositionsBackup.length; e++) this.vertices[e].position = this.vertexPositionsBackup[e], this.vertices[e].positioned = !0; for (var t = 0; t < this.ringPositionsBackup.length; t++) this.rings[t].center = this.ringPositionsBackup[t], this.rings[t].positioned = !0 } }, { key: "backupRingInformation", value: function() { this.originalRings = [], this.originalRingConnections = []; for (var e = 0; e < this.rings.length; e++) this.originalRings.push(this.rings[e]); for (var t = 0; t < this.ringConnections.length; t++) this.originalRingConnections.push(this.ringConnections[t]); for (var i = 0; i < this.vertices.length; i++) this.vertices[i].value.backupRings() } }, { key: "restoreRingInformation", value: function() { var e = this.getBridgedRings();
                this.rings = [], this.ringConnections = []; for (var t = 0; t < e.length; t++)
                    for (var i = e[t], r = 0; r < i.rings.length; r++) { var n = i.rings[r];
                        this.originalRings[n.id].center = n.center }
                for (var s = 0; s < this.originalRings.length; s++) this.rings.push(this.originalRings[s]); for (var o = 0; o < this.originalRingConnections.length; o++) this.ringConnections.push(this.originalRingConnections[o]); for (var a = 0; a < this.vertices.length; a++) this.vertices[a].value.restoreRings() } }, { key: "createRing", value: function(e) { var t = arguments.length > 1 && void 0 !== arguments[1] ? arguments[1] : null,
                    i = arguments.length > 2 && void 0 !== arguments[2] ? arguments[2] : null,
                    r = this,
                    n = arguments.length > 3 && void 0 !== arguments[3] ? arguments[3] : null,
                    s = arguments.length > 4 && void 0 !== arguments[4] && arguments[4]; if (!e.positioned || s) { t = t || new Vector2(0, 0); var o = e.getOrderedNeighbours(this.ringConnections),
                        a = i ? Vector2.subtract(i.position, t).angle() : 0,
                        l = MathHelper.polyCircumradius(this.opts.bondLength, e.getSize()),
                        h = MathHelper.centralAngle(e.getSize());
                    e.centralAngle = h; var u = a,
                        g = this; if (!s) { if (e.eachMember(this.vertices, function(i) { var r = g.vertices[i];
                                r.positioned || (r.position.x = t.x + Math.cos(u) * l, r.position.y = t.y + Math.sin(u) * l), u += h, (!e.isBridged || e.rings.length < 3) && (r.positioned = !0) }, i ? i.id : null, n ? n.id : null), e.isBridged) { var c = ArrayHelper.merge(e.members, e.insiders);
                            this.forceLayout(c, t, i.id, e) }
                        this.vertices[e.members[0]].value.addAnchoredRing(e.id), e.positioned = !0, e.center = t } for (var v = 0; v < o.length; v++) { var d = this.getRing(o[v].neighbour); if (!d.positioned) { var f = RingConnection.getVertices(this.ringConnections, e.id, d.id); if (2 == f.length) ! function() { e.isFused = !0, d.isFused = !0; var t = r.vertices[f[0]],
                                    i = r.vertices[f[1]],
                                    n = Vector2.midpoint(t.position, i.position),
                                    s = Vector2.normals(t.position, i.position);
                                ArrayHelper.each(s, function(e) { e.normalize() }); var o = MathHelper.polyCircumradius(r.opts.bondLength, d.getSize()),
                                    a = MathHelper.apothem(o, d.getSize());
                                ArrayHelper.each(s, function(e) { e.multiply(a) }), ArrayHelper.each(s, function(e) { e.add(n) }); var l = s[0];
                                r.isPointInRing(l) && (l = s[1]); var h = Vector2.subtract(t.position, l),
                                    u = Vector2.subtract(i.position, l); - 1 === h.clockwise(u) ? r.createRing(d, l, t, i) : r.createRing(d, l, i, t) }();
                            else if (1 == f.length) { e.isSpiro = !0, d.isSpiro = !0; var p = this.vertices[f[0]],
                                    y = Vector2.subtract(t, p.position);
                                y.invert(), y.normalize(); var m = MathHelper.polyCircumradius(this.opts.bondLength, d.getSize());
                                y.multiply(m), y.add(p.position), this.createRing(d, y, p) } } } for (var b = 0; b < e.members.length; b++)
                        for (var A = this.vertices[e.members[b]], x = A.getNeighbours(), k = 0; k < x.length; k++)
                            if (!e.thisOrNeighboursContain(this.rings, x[k])) { var C = this.vertices[x[k]];
                                this.createNextBond(C, A, e.center) } } } }, { key: "rotateSubtree", value: function(e, t, i, r) { var n = this;
                this.traverseTree(e, t, function(e) { e.position.rotateAround(i, r); for (var t = 0; t < e.value.anchoredRings.length; t++) { var s = n.rings[e.value.anchoredRings[t]];
                        s && s.center.rotateAround(i, r) } }) } }, { key: "getSubtreeOverlapScore", value: function(e, t, i) { var r = this,
                    n = 0,
                    s = new Vector2; return this.traverseTree(e, t, function(e) { var t = i[e.id];
                    n += t; var o = r.vertices[e.id].position.clone();
                    o.multiply(t), s.add(o) }), s.divide(n), { value: n, center: s } } }, { key: "getCurrentCenterOfMass", value: function() { for (var e = new Vector2, t = 0, i = 0; i < this.vertices.length; i++) { var r = this.vertices[i];
                    r.positioned && (e.add(r.position), t++) } return e.divide(t) } }, {
            key: "getCurrentCenterOfMassInNeigbourhood",
            value: function(e) {
                for (var t = arguments.length > 1 && void 0 !== arguments[1] ? arguments[1] : 2 * this.opts.bondLength, i = new Vector2, r = 0, n = t * t, s = 0; s < this.vertices.length; s++) { var o = this.vertices[s];
                    o.positioned && e.distanceSq(o.position) < n && (i.add(o.position), r++) }
                return i.divide(r)
            }
        }, { key: "resolvePrimaryOverlaps", value: function() { for (var e = [], t = [], i = new Array(this.vertices.length), r = 0; r < this.rings.length; r++)
                    for (var n = this.rings[r], s = 0; s < n.members.length; s++) { var o = this.vertices[n.members[s]]; if (!i[o.id] && (i[o.id] = !0, o.getNeighbours().length > 2)) { for (var a = [], l = 0; l < o.value.rings.length; l++) a.push(o.value.rings[l]);
                            e.push({ common: o, rings: a, vertices: this.getNonRingNeighbours(o.id) }) } }
                for (var h = 0; h < t.length; h++) { var u = t[h],
                        g = -u.vertex.position.getRotateToAngle(u.other.position, u.common.position);
                    this.rotateSubtree(u.vertex.id, u.common.id, g + Math.PI, u.common.position) } for (var c = 0; c < e.length; c++) { var v = e[c]; if (1 == v.vertices.length) { var d = v.vertices[0]; if (1 == d.getNeighbours().length) { d.flippable = !0, d.flipCenter = v.common.id; for (var f = 0; f < v.rings.length; f++) d.flipRings.push(v.rings[f]) } if (2 === v.rings.length) { for (var p = v.common.getNeighbours(), y = [], m = 0; m < p.length; m++) { var b = this.vertices[p[m]];
                                this.isRingConnection(b.id, v.common.id) || b.id === d.id || y.push(b.position) } var A = Vector2.midpoint(y[0], y[1]),
                                x = d.position.getRotateToAngle(A, v.common.position);
                            x *= d.position.clockwise(A), this.rotateSubtree(d.id, v.common.id, x, v.common.position) } } else if (2 == v.vertices.length) { var k = (2 * Math.PI - this.getRing(v.rings[0]).getAngle()) / 6,
                            C = v.vertices[0],
                            R = v.vertices[1]; if (C.backAngle -= k, R.backAngle += k, this.rotateSubtree(C.id, v.common.id, k, v.common.position), this.rotateSubtree(R.id, v.common.id, -k, v.common.position), 1 == C.getNeighbours().length) { C.flippable = !0, C.flipCenter = v.common.id, C.flipNeighbour = R.id; for (var w = 0; w < v.rings.length; w++) C.flipRings.push(v.rings[w]) } if (1 == R.getNeighbours().length) { R.flippable = !0, R.flipCenter = v.common.id, R.flipNeighbour = C.id; for (var S = 0; S < v.rings.length; S++) R.flipRings.push(v.rings[S]) } } } } }, { key: "resolveSecondaryOverlaps", value: function(e) { for (var t = 0; t < e.length; t++)
                    if (e[t].score > this.opts.bondLength / (4 * this.opts.bondLength)) { var i = this.vertices[e[t].id]; if (i.isTerminal()) { var r = this.getClosestVertex(i); if (r) { var n = null;
                                n = r.isTerminal() ? 0 === r.id ? this.vertices[1].position : r.previousPosition : 0 === r.id ? this.vertices[1].position : r.position; var s = 0 === i.id ? this.vertices[1].position : i.previousPosition;
                                i.position.rotateAwayFrom(n, s, MathHelper.toRad(20)) } } if (i.flippable) { var o = i.flipRings[0] ? this.rings[i.flipRings[0]] : null,
                                a = i.flipRings[1] ? this.rings[i.flipRings[1]] : null,
                                l = this.vertices[i.flipCenter].position; if (o && a) { var h = o.members.length > a.members.length ? o : a;
                                a = o.members.length < a.members.length ? o : a, o = h }
                            this.opts.allowFlips && (o && o.allowsFlip() ? (i.position.rotateTo(o.center, l), o.setFlipped(), i.flipNeighbour) : a && a.allowsFlip() && (i.position.rotateTo(a.center, l), a.setFlipped(), i.flipNeighbour)) } } } }, { key: "createNextBond", value: function(e, t, i, r) { if (!e.positioned) { if (t)
                        if (0 !== t.value.rings.length || e.value.isBridge || t.value.isBridge) { if (t.value.isBridgeNode && e.value.isBridge) { var n = Vector2.subtract(i, t.position);
                                n.normalize(), n.multiply(this.opts.bondLength), e.position.add(t.position), e.position.add(n), e.previousPosition = t.position, e.positioned = !0 } else if (e.value.isBridge) { var s = new Vector2(this.opts.bondLength, 0);
                                s.rotate(i), s.add(t.position), e.globalAngle = i, e.position = s, e.previousPosition = t.position, e.positioned = !0 } else if (1 === t.value.rings.length || t.value.isBridge) { var o = Vector2.subtract(i, t.position);
                                o.invert(), o.normalize(), o.multiply(this.opts.bondLength), e.position.add(t.position), e.position.add(o), e.previousPosition = t.position, e.positioned = !0 } else if (2 == t.value.rings.length) { var a = this.getRing(t.value.rings[0]),
                                    l = this.getRing(t.value.rings[1]),
                                    h = Vector2.subtract(l.center, a.center),
                                    u = Vector2.subtract(t.position, a.center),
                                    g = Vector2.scalarProjection(u, h);
                                h.normalize(), h.multiply(g), h.add(a.center); var c = Vector2.subtract(h, t.position);
                                c.invert(), c.normalize(), c.multiply(this.opts.bondLength), e.position.add(t.position), e.position.add(c), e.previousPosition = t.position, e.positioned = !0 } } else { var v = new Vector2(this.opts.bondLength, 0);
                            v.rotate(i), v.add(t.position), e.globalAngle = i, e.position = v, e.previousPosition = t.position, e.positioned = !0 }
                    else { var d = new Vector2(this.opts.bondLength, 0);
                        d.rotate(MathHelper.toRad(-120)), e.previousPosition = d, e.position = new Vector2(this.opts.bondLength, 0), e.angle = MathHelper.toRad(-120), e.globalAngle = e.angle, e.positioned = !0 } if (e.value.rings.length > 0) { var f = this.getRing(e.value.rings[0]),
                            p = Vector2.subtract(e.previousPosition, e.position);
                        p.invert(), p.normalize(); var y = MathHelper.polyCircumradius(this.opts.bondLength, f.getSize());
                        p.multiply(y), p.add(e.position), this.createRing(f, p, e) } else { var m = e.getNeighbours();
                        t && (m = ArrayHelper.remove(m, t.id)); var b = e.getAngle(); if (1 === m.length) { var A = this.vertices[m[0]]; if ("#" === e.value.bondType || t && "#" === t.value.bondType || "=" === e.value.bondType && t && "=" === t.value.bondType) { if (e.value.explicit = !0, t) { this.getEdge(e.id, t.id).center = !0 }
                                this.getEdge(e.id, A.id).center = !0, A.globalAngle = b, A.angle = 0, this.createNextBond(A, e, A.globalAngle, -r) } else if (t && t.value.rings.length > 0) { var x = MathHelper.toRad(60),
                                    k = -x,
                                    C = new Vector2(this.opts.bondLength, 0),
                                    R = new Vector2(this.opts.bondLength, 0);
                                C.rotate(x).add(e.position), R.rotate(k).add(e.position); var w = this.getCurrentCenterOfMass(),
                                    S = C.distance(w),
                                    V = R.distance(w);
                                A.angle = S < V ? k : x, r = A.angle > 0 ? -1 : 1, A.globalAngle = b + A.angle, this.createNextBond(A, e, A.globalAngle, r) } else { if (r) A.angle = MathHelper.toRad(60) * r, r = -r;
                                else { var N = MathHelper.toRad(60),
                                        L = -N,
                                        B = new Vector2(this.opts.bondLength, 0),
                                        T = new Vector2(this.opts.bondLength, 0);
                                    B.rotate(N).add(e.position), T.rotate(L).add(e.position); var I = this.getCurrentCenterOfMass(),
                                        H = B.distance(I),
                                        M = T.distance(I);
                                    A.angle = H < M ? L : N, r = A.angle > 0 ? -1 : 1 }
                                A.globalAngle = b + A.angle, this.createNextBond(A, e, A.globalAngle, r) } } else if (2 === m.length) { var P = this.getTreeDepth(m[0], e.id),
                                O = this.getTreeDepth(m[1], e.id),
                                W = 0,
                                F = 1; if (P > O && (W = 1, F = 0), 1 === e.position.clockwise(e.previousPosition)) { var E = this.vertices[m[W]],
                                    D = this.vertices[m[F]];
                                D.angle = MathHelper.toRad(60), E.angle = -MathHelper.toRad(60), D.globalAngle = b + D.angle, E.globalAngle = b + E.angle, this.createNextBond(D, e, D.globalAngle, -r), this.createNextBond(E, e, E.globalAngle, -r) } else { var z = this.vertices[m[W]],
                                    _ = this.vertices[m[F]];
                                _.angle = -MathHelper.toRad(60), z.angle = MathHelper.toRad(60), _.globalAngle = b + _.angle, z.globalAngle = b + z.angle, this.createNextBond(z, e, z.globalAngle, -r), this.createNextBond(_, e, _.globalAngle, -r) } } else if (3 === m.length) { var q = this.getTreeDepth(m[0], e.id),
                                U = this.getTreeDepth(m[1], e.id),
                                j = this.getTreeDepth(m[2], e.id),
                                X = this.vertices[m[0]],
                                G = this.vertices[m[1]],
                                Y = this.vertices[m[2]]; if (U > q && U > j ? (X = this.vertices[m[1]], G = this.vertices[m[0]], Y = this.vertices[m[2]]) : j > q && j > U && (X = this.vertices[m[2]], G = this.vertices[m[0]], Y = this.vertices[m[1]]), 1 === this.getTreeDepth(G.id, e.id) && 1 === this.getTreeDepth(Y.id, e.id) && this.getTreeDepth(X.id, e.id) > 1) { if (r) X.angle = MathHelper.toRad(60) * r, r = -r;
                                else { var Z = MathHelper.toRad(60),
                                        K = -Z,
                                        $ = new Vector2(this.opts.bondLength, 0),
                                        J = new Vector2(this.opts.bondLength, 0);
                                    $.rotate(Z).add(e.position), J.rotate(K).add(e.position); var Q = this.getCurrentCenterOfMass(),
                                        ee = $.distance(Q),
                                        te = J.distance(Q);
                                    X.angle = ee < te ? K : Z, r = X.angle > 0 ? -1 : 1 }
                                X.globalAngle = b + X.angle, this.createNextBond(X, e, X.globalAngle, -r), e.value.bracket && "@@" === e.value.bracket.chirality ? (Y.angle = MathHelper.toRad(30) * r, G.angle = MathHelper.toRad(90) * r, Y.globalAngle = b + Y.angle, G.globalAngle = b + G.angle, this.createNextBond(Y, e, Y.globalAngle), this.createNextBond(G, e, G.globalAngle)) : (G.angle = MathHelper.toRad(30) * r, Y.angle = MathHelper.toRad(90) * r, G.globalAngle = b + G.angle, Y.globalAngle = b + Y.angle, this.createNextBond(G, e, G.globalAngle), this.createNextBond(Y, e, Y.globalAngle)) } else X.angle = 0, G.angle = MathHelper.toRad(90), Y.angle = -MathHelper.toRad(90), X.globalAngle = b + X.angle, G.globalAngle = b + G.angle, Y.globalAngle = b + Y.angle, this.createNextBond(X, e, X.globalAngle), this.createNextBond(G, e, G.globalAngle), this.createNextBond(Y, e, Y.globalAngle) } else if (4 === m.length) { var ie = this.getTreeDepth(m[0], e.id),
                                re = this.getTreeDepth(m[1], e.id),
                                ne = this.getTreeDepth(m[2], e.id),
                                se = this.getTreeDepth(m[3], e.id),
                                oe = this.vertices[m[0]],
                                ae = this.vertices[m[1]],
                                le = this.vertices[m[2]],
                                he = this.vertices[m[3]];
                            re > ie && re > ne && re > se ? (oe = this.vertices[m[1]], ae = this.vertices[m[0]], le = this.vertices[m[2]], he = this.vertices[m[3]]) : ne > ie && ne > re && ne > se ? (oe = this.vertices[m[2]], ae = this.vertices[m[0]], le = this.vertices[m[1]], he = this.vertices[m[3]]) : se > ie && se > re && se > ne && (oe = this.vertices[m[3]], ae = this.vertices[m[0]], le = this.vertices[m[1]], he = this.vertices[m[2]]), oe.angle = -MathHelper.toRad(36), ae.angle = MathHelper.toRad(36), le.angle = -MathHelper.toRad(108), he.angle = MathHelper.toRad(108), oe.globalAngle = b + oe.angle, ae.globalAngle = b + ae.angle, le.globalAngle = b + le.angle, he.globalAngle = b + he.angle, this.createNextBond(oe, e, oe.globalAngle), this.createNextBond(ae, e, ae.globalAngle), this.createNextBond(le, e, le.globalAngle), this.createNextBond(he, e, he.globalAngle) } } } } }, { key: "getCommonRingbondNeighbour", value: function(e) { for (var t = e.getNeighbours(), i = 0; i < t.length; i++) { var r = this.vertices[t[i]]; if (ArrayHelper.containsAll(r.value.rings, e.value.rings)) return r } return null } }, { key: "isPointInRing", value: function(e) { for (var t = 0; t < this.rings.length; t++) { var i = this.rings[t]; if (i.positioned) { var r = (i.getPolygon(this.vertices), MathHelper.polyCircumradius(this.opts.bondLength, i.getSize())),
                            n = r * r; if (e.distanceSq(i.center) < n) return !0 } } return !1 } }, { key: "isEdgeInRing", value: function(e) { var t = this.vertices[e.sourceId],
                    i = this.vertices[e.targetId]; return this.areVerticesInSameRing(t, i) } }, { key: "isEdgeRotatable", value: function(e) { var t = this.vertices[e.sourceId],
                    i = this.vertices[e.targetId]; return "-" === e.bondType && (!(t.getNeighbourCount() + i.getNeighbourCount() < 5) && !(t.value.rings.length > 0 && i.value.rings.length > 0 && this.areVerticesInSameRing(t, i))) } }, { key: "isRingAromatic", value: function(e) { for (var t = 0; t < e.members.length; t++)
                    if (!this.isVertexInAromaticRing(e.members[t])) return !1;
                return !0 } }, { key: "isEdgeInAromaticRing", value: function(e) { return this.isVertexInAromaticRing(e.sourceId) && this.isVertexInAromaticRing(e.targetId) } }, { key: "isVertexInAromaticRing", value: function(e) { var t = this.vertices[e].value.element; return t == t.toLowerCase() } }, { key: "getEdgeNormals", value: function(e) { var t = this.vertices[e.sourceId].position,
                    i = this.vertices[e.targetId].position,
                    r = Vector2.normals(t, i); return ArrayHelper.each(r, function(e) { e.normalize() }), r } }, { key: "getTreeDepth", value: function(e, t) { for (var i = this.vertices[e].getSpanningTreeNeighbours(t), r = 0, n = 0; n < i.length; n++) { var s = i[n],
                        o = this.getTreeDepth(s, e);
                    o > r && (r = o) } return r + 1 } }, { key: "traverseTree", value: function(e, t, i) { var r = arguments.length > 3 && void 0 !== arguments[3] ? arguments[3] : null,
                    n = arguments.length > 4 && void 0 !== arguments[4] && arguments[4],
                    s = arguments.length > 5 && void 0 !== arguments[5] ? arguments[5] : 1,
                    o = arguments.length > 6 && void 0 !== arguments[6] ? arguments[6] : []; if (!(null !== r && s > r + 1)) { for (var a = 0; a < o.length; a++)
                        if (o[a] === e) return;
                    o.push(e); var l = this.vertices[e],
                        h = l.getNeighbours(t);
                    (!n || s > 1) && i(l); for (var u = 0; u < h.length; u++) this.traverseTree(h[u], e, i, r, n, s + 1, o) } } }, { key: "getBondCount", value: function(e) { for (var t = 0, i = 0; i < e.edges.length; i++) t += this.edges[e.edges[i]].getBondCount(); return t } }, { key: "getNonRingNeighbours", value: function(e) { for (var t = [], i = this.vertices[e], r = i.getNeighbours(), n = 0; n < r.length; n++) { var s = this.vertices[r[n]];
                    0 === ArrayHelper.intersection(i.value.rings, s.value.rings).length && 0 == s.value.isBridge && t.push(s) } return t } }, { key: "annotateChirality", value: function() { for (var e = 0; e < this.vertices.length; e++) { var t = this.vertices[e]; if (t.value.bracket && "c" === t.value.element.toLowerCase() && 4 === t.getNeighbours().length || t.value.bracket && t.value.bracket.hcount > 0 && 3 === t.getNeighbours().length) { var i = t.value.bracket.chirality; if (null === i) continue; var r = t.getNeighbours(),
                            n = new Array(r.length);
                        this.vertices[t.parentVertexId].value.setOrder(t.id, 0); for (var s = 0; s < r.length; s++) { var o = r[s];
                            o !== t.parentVertexId ? n[this.vertices[o].value.getOrder(t.id)] = o : n[0] = o } if ("@" === i) { var a = this.getEdge(n[3], t.id); "down" !== a.chiral && (a.chiral = "up"); var l = this.getEdge(n[1], t.id); "up" !== l.chiral && (l.chiral = "down") } else if ("@@" === i) { var h = this.getEdge(n[1], t.id); "down" !== h.chiral && (h.chiral = "up"); var u = this.getEdge(n[3], t.id); "up" !== u.chiral && (u.chiral = "down") } } } } }], [{ key: "clean", value: function(e) { return e.replace(/[^A-Za-z0-9@\.\+\-\?!\(\)\[\]\{\}\/\\=#\$:\*]/g, "") } }, { key: "apply", value: function(t) { for (var i = arguments.length > 1 && void 0 !== arguments[1] ? arguments[1] : "light", r = new e(t), n = document.querySelectorAll("canvas[data-smiles]"), s = 0; s < n.length; s++) { var o = n[s],
                        a = e.parse(e.clean(o.getAttribute("data-smiles")));
                    r.draw(a, o, i, !1) } } }, { key: "parse", value: function(e) { return SMILESPARSER.parse(e) } }]), e
    }(),
    Vector2 = function() {
        function e(t, i) { _classCallCheck(this, e), 0 == arguments.length ? (this.x = 0, this.y = 0) : 1 == arguments.length ? (this.x = t.x, this.y = t.y) : (this.x = t, this.y = i) } return _createClass(e, [{ key: "set", value: function() { var e = arguments.length > 0 && void 0 !== arguments[0] ? arguments[0] : 0,
                    t = arguments.length > 1 && void 0 !== arguments[1] ? arguments[1] : 0; return this.x = e, this.y = t, this } }, { key: "clone", value: function() { return new e(this.x, this.y) } }, { key: "toString", value: function() { return "(" + this.x + "," + this.y + ")" } }, { key: "add", value: function(e) { return this.x += e.x, this.y += e.y, this } }, { key: "subtract", value: function(e) { return this.x -= e.x, this.y -= e.y, this } }, { key: "divide", value: function(e) { return this.x /= e, this.y /= e, this } }, { key: "multiply", value: function(e) { return this.x *= e, this.y *= e, this } }, { key: "invert", value: function() { return this.x = -this.x, this.y = -this.y, this } }, { key: "angle", value: function() { return Math.atan2(this.y, this.x) } }, { key: "distance", value: function(e) { return Math.sqrt((e.x - this.x) * (e.x - this.x) + (e.y - this.y) * (e.y - this.y)) } }, { key: "distanceSq", value: function(e) { return (e.x - this.x) * (e.x - this.x) + (e.y - this.y) * (e.y - this.y) } }, { key: "clockwise", value: function(e) { var t = this.y * e.x,
                    i = this.x * e.y; return t > i ? -1 : t === i ? 0 : 1 } }, { key: "rotate", value: function(t) { var i = new e; return i.x = this.x * Math.cos(t) - this.y * Math.sin(t), i.y = this.x * Math.sin(t) + this.y * Math.cos(t), this.x = i.x, this.y = i.y, this } }, { key: "rotateAround", value: function(e, t) { var i = Math.sin(e),
                    r = Math.cos(e);
                this.x -= t.x, this.y -= t.y; var n = this.x * r - this.y * i,
                    s = this.x * i + this.y * r; return this.x = n + t.x, this.y = s + t.y, this } }, { key: "rotateTo", value: function(t, i) { var r = arguments.length > 2 && void 0 !== arguments[2] ? arguments[2] : 0;
                this.x += .001, this.y -= .001; var n = e.subtract(this, i),
                    s = e.subtract(t, i),
                    o = e.angle(s, n); return this.rotateAround(o + r, i), this } }, { key: "rotateAwayFrom", value: function(e, t, i) { this.rotateAround(i, t); var r = this.distanceSq(e);
                this.rotateAround(-2 * i, t), this.distanceSq(e) < r && this.rotateAround(2 * i, t) } }, { key: "getRotateAwayFromAngle", value: function(e, t, i) { var r = this.clone();
                r.rotateAround(i, t); var n = r.distanceSq(e); return r.rotateAround(-2 * i, t), r.distanceSq(e) < n ? i : -i } }, { key: "getRotateTowardsAngle", value: function(e, t, i) { var r = this.clone();
                r.rotateAround(i, t); var n = r.distanceSq(e); return r.rotateAround(-2 * i, t), r.distanceSq(e) > n ? i : -i } }, { key: "getRotateToAngle", value: function(t, i) { var r = e.subtract(this, i),
                    n = e.subtract(t, i),
                    s = e.angle(n, r); return Number.isNaN(s) ? 0 : s } }, { key: "isInPolygon", value: function(e) { for (var t = !1, i = 0, r = e.length - 1; i < e.length; r = i++) e[i].y > this.y != e[r].y > this.y && this.x < (e[r].x - e[i].x) * (this.y - e[i].y) / (e[r].y - e[i].y) + e[i].x && (t = !t); return t } }, { key: "length", value: function() { return Math.sqrt(this.x * this.x + this.y * this.y) } }, { key: "normalize", value: function() { return this.divide(this.length()), this } }, { key: "normalized", value: function() { return e.divide(this, this.length()) } }, { key: "whichSide", value: function(e, t) { return (this.x - e.x) * (t.y - e.y) - (this.y - e.y) * (t.x - e.x) } }, { key: "sameSideAs", value: function(e, t, i) { var r = this.whichSide(e, t),
                    n = i.whichSide(e, t); return r < 0 && n < 0 || 0 == r && 0 == n || r > 0 && n > 0 } }], [{ key: "add", value: function(t, i) { return new e(t.x + i.x, t.y + i.y) } }, { key: "subtract", value: function(t, i) { return new e(t.x - i.x, t.y - i.y) } }, { key: "multiply", value: function(t, i) { return i.x && i.y ? new e(t.x * i.x, t.y * i.y) : new e(t.x * i, t.y * i) } }, { key: "multiplyScalar", value: function(t, i) { return new e(t).multiply(i) } }, { key: "midpoint", value: function(t, i) { return new e((t.x + i.x) / 2, (t.y + i.y) / 2) } }, { key: "normals", value: function(t, i) { var r = e.subtract(i, t); return [new e(-r.y, r.x), new e(r.y, -r.x)] } }, { key: "divide", value: function(t, i) { return i.x && i.y ? new e(t.x / i.x, t.y / i.y) : new e(t.x / i, t.y / i) } }, { key: "dot", value: function(e, t) { return e.x * t.x + e.y * t.y } }, { key: "angle", value: function(t, i) { var r = e.dot(t, i); return Math.acos(r / (t.length() * i.length())) } }, { key: "threePointangle", value: function(t, i, r) { var n = e.subtract(i - t),
                    s = e.subtract(r, i),
                    o = (t.distance(i), i.distance(r)); return Math.acos(e.dot(n, s) / (abLenght * o)) } }, { key: "scalarProjection", value: function(t, i) { var r = i.normalized(); return e.dot(t, r) } }]), e }(),
    Vertex = function() {
        function e(t) { var i = arguments.length > 1 && void 0 !== arguments[1] ? arguments[1] : 0,
                r = arguments.length > 2 && void 0 !== arguments[2] ? arguments[2] : 0;
            _classCallCheck(this, e), this.id = null, this.value = t, this.position = new Vector2(i || 0, r || 0), this.previousPosition = new Vector2(0, 0), this.parentVertexId = null, this.children = [], this.spanningTreeChildren = [], this.edges = [], this.positioned = !1, this.angle = 0, this.globalAngle = 0, this.dir = 1, this.backAngle = 0, this.flippable = !1, this.flipCenter = null, this.flipNeighbour = null, this.flipRings = new Array } return _createClass(e, [{ key: "isTerminal", value: function() { return null === this.parentVertexId && this.children.length < 2 || 0 === this.children.length } }, { key: "clone", value: function() { var t = new e(this.value, this.position.x, this.position.y); return t.id = this.id, t.previousPosition = new Vector2(this.previousPosition.x, this.previousPosition.y), t.parentVertexId = this.parentVertexId, t.children = ArrayHelper.clone(this.children), t.spanningTreeChildren = ArrayHelper.clone(this.spanningTreeChildren), t.edges = ArrayHelper.clone(this.edges), t.positioned = this.positioned, t.angle = this.angle, t.backAngle = this.backAngle, t.flippable = this.flippable, t.flipCenter = this.flipCenter, t.flipRings = ArrayHelper.clone(this.flipRings), t } }, { key: "equals", value: function(e) { return this.id === e.id } }, { key: "getAngle", value: function() { var e = arguments.length > 0 && void 0 !== arguments[0] ? arguments[0] : null,
                    t = arguments.length > 1 && void 0 !== arguments[1] && arguments[1],
                    i = null; return i = e ? Vector2.subtract(this.position, e) : Vector2.subtract(this.position, this.previousPosition), t ? MathHelper.toDeg(i.angle()) : i.angle() } }, { key: "getTextDirection", value: function(e) { for (var t = this.getNeighbours(), i = [], r = 0; r < t.length; r++) i.push(this.getAngle(e[t[r]].position)); var n = MathHelper.meanAngle(i),
                    s = Math.PI / 2; return n = Math.round(Math.round(n / s) * s, 3), 2 == n ? "down" : -2 == n ? "up" : 0 === n || -0 === n ? "right" : 3 == n || -3 == n ? "left" : "down" } }, { key: "getNeighbours", value: function() { for (var e = arguments.length > 0 && void 0 !== arguments[0] ? arguments[0] : null, t = [], i = 0; i < this.children.length; i++) void 0 !== e && e == this.children[i] || t.push(this.children[i]); return null != this.parentVertexId && (void 0 !== e && e == this.parentVertexId || t.push(this.parentVertexId)), t } }, { key: "getNeighbourCount", value: function() { var e = this.children.length; return null !== this.parentVertexId && (e += 1), e } }, { key: "getCommonNeighbours", value: function(e) { for (var t = new Array, i = this.getNeighbours(), r = e.getNeighbours(), n = 0; n < i.length; n++)
                    for (var s = 0; s < r.length; s++) i[n] === r[s] && t.push(i[n]); return t } }, { key: "isNeighbour", value: function(e) { if (this.parentVertexId === e) return !0; for (var t = 0; t < this.children.length; t++)
                    if (this.children[t] === e) return !0 } }, { key: "getSpanningTreeNeighbours", value: function() { for (var e = arguments.length > 0 && void 0 !== arguments[0] ? arguments[0] : null, t = [], i = 0; i < this.spanningTreeChildren.length; i++) void 0 !== e && e == this.spanningTreeChildren[i] || t.push(this.spanningTreeChildren[i]); return null != this.parentVertexId && (void 0 !== e && e == this.parentVertexId || t.push(this.parentVertexId)), t } }, { key: "getNextInRing", value: function(e, t, i) { for (var r = this.getNeighbours(), n = 0; n < r.length; n++)
                    if (ArrayHelper.contains(e[r[n]].value.rings, { value: t }) && r[n] != i) return r[n];
                return null } }]), e }();